#' Simulate data from the categorical spatial mark resight model
#' @description Will document after I defend. Should be able to get the idea in the example found in the mcmc.conCatSMR.ID help file.
#' @author Ben Augustine
#' @export
sim.conCatSMR.ID <-
  function(N=50,n.marked=10,lam0=0.25,sigma=0.50,K=10,X=X,buff=3,obstype="bernoulli",ncat=ncat,
           pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=c(1,1),tlocs=0,pID=1,
           marktype="premarked"){
    if(!marktype%in%c("natural","premarked")){
      stop("marktype must be 'natural' or 'premarked'")
    }
    library(abind)
    # simulate a population of activity centers
    X=as.matrix(X)
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J=nrow(X)
    #simulate IDcovs
    G.true=matrix(NA,nrow=N,ncol=ncat) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:ncat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    
    # Capture and mark individuals
    y.sight <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.sight[i,j,k]=rbinom(1,1,pd[i,j]) 
          }
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.sight[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    if(marktype=="natural"){
      #reorder data so enough marked guys are at the top
      #must be random sample of marked guys, e.g. not ordered by # of captures
      idx=which(rowSums(y.sight)>0)
      move=sample(idx,n.marked)
      idx=setdiff(1:N,move)
      y.sight=y.sight[c(move,idx),,]
      s=s[c(move,idx),]
    }
    IDmarked=1:n.marked
    umguys=setdiff(1:N,IDmarked)
    
    #split sightings into marked and unmarked histories, considering occasion of marking
    # y.sightMK=apply(y.sight,c(1,3),sum)
    # n.marked=sum(rowSums(markedS)!=0)
    y.sight.marked=y.sight[IDmarked,,]
    # y.sight.unmarked=y.sight[IDum,,]
    G.marked=G.true[IDmarked,]
    n.samples=sum(y.sight[umguys,,])
    G.unmarked=matrix(NA,nrow=n.samples,ncol=ncat)
    y.sight.unmarked=array(0,dim=c(n.samples,J,K))
    IDum=rep(NA,n.samples)
    idx=1
    for(i in 1:length(umguys)){
      for(j in 1:J){ #then traps
        for(k in 1:K){ #then occasions
          if(y.sight[umguys[i],j,k]>0){ #is there at least one sample here?
              for(l in 1:y.sight[umguys[i],j,k]){ #then samples
                G.unmarked[idx,]=G.true[i,]
                y.sight.unmarked[idx,j,k]=1
                IDum[idx]=umguys[i]
                idx=idx+1
              }
            }
          }
        }
    }
    if(!is.matrix(G.marked)){
      G.marked=matrix(G.marked,ncol=1)
    }
    if(!is.matrix(G.unmarked)){
      G.unmarked=matrix(G.unmarked,ncol=1)
    }
    #Mark status ID process
    if(pMarkID[2]<1&pMarkID[1]==1){#can perfectly ID marked status of marked guys but not unmarked
      unm2unk=which(rbinom(n.samples,1,pMarkID[2])==0)
      if(length(unm2unk)>0){
        G.unk=G.unmarked[unm2unk,]
        G.unmarked=G.unmarked[-unm2unk,]
        if(!is.matrix(G.unk)){
          G.unk=matrix(G.unk,ncol=1)
        }
        if(!is.matrix(G.unmarked)){
          G.unmarked=matrix(G.unmarked,ncol=1)
        }
        y.sight.unk=y.sight.unmarked[unm2unk,,]
        y.sight.unmarked=y.sight.unmarked[-unm2unk,,]
        IDunk=IDum[unm2unk]
        IDum=IDum[-unm2unk]
      }else{
        y.sight.unk=NA
        G.unk=NA
        IDunk=NA
      }
    }else if(pMarkID[2]<1&pMarkID[1]<1){#Can't perfectly ID mark status of marked or unmarked guys
      #unmarked guys
      unm2unk=which(rbinom(n.samples,1,pMarkID[2])==0)
      if(length(unm2unk)>0){
        G.unk=G.unmarked[unm2unk,]
        G.unmarked=G.unmarked[-unm2unk,]
        if(!is.matrix(G.unk)){
          G.unk=matrix(G.unk,ncol=1)
        }
        if(!is.matrix(G.unmarked)){
          G.unmarked=matrix(G.unmarked,ncol=1)
        }
        y.sight.unk=y.sight.unmarked[unm2unk,,]
        y.sight.unmarked=y.sight.unmarked[-unm2unk,,]
        IDunk=IDum[unm2unk]
        IDum=IDum[-unm2unk]
      }else{#didn't lose any mark status info
        G.unk=matrix(0,nrow=0,ncol=ncat)
        y.sight.unmarked=array(0,dim=c(0,J,K))
        IDunk=rep(0,0)
      }
      #add marked guys
      idx1=which(y.sight.marked>0)#used to extract counts
      count=y.sight.marked[idx1]
      idx2=which(y.sight.marked>0,arr.ind=TRUE)#used to move counts
      idx3=rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2=idx2[idx3,]
      rem=1-rbinom(nrow(idx2),1,pMarkID[1])#which counts can't we tell marked status
      if(sum(rem)>0){
        idx2=idx2[which(rem==1),]
        if(!is.matrix(idx2)){
          idx2=matrix(idx2,nrow=1)
        }
        y.sight.unk2=array(0,dim=c(nrow(idx2),J,K))
        G.unk2=matrix(NA,nrow=nrow(idx2),ncol=ncat)
        for(i in 1:nrow(idx2)){
          #remove unk guys from marked sightings
          y.sight.marked[idx2[i,1],idx2[i,2],idx2[i,3]]=y.sight.marked[idx2[i,1],idx2[i,2],idx2[i,3]]-1
          #add unk guys to new unk structure
          y.sight.unk2[i,idx2[i,2],idx2[i,3]]=1
          #add unk guy IDcovs to new structure
          G.unk2[i,]=G.marked[idx2[i,1],]
        }
        #paste new unk structures to old ones
        y.sight.unk=abind(y.sight.unk,y.sight.unk2,along=1)
        G.unk=rbind(G.unk,G.unk2)
        IDunk=c(IDunk,idx2[,1])
      }
      if(dim(G.unk)[1]==0){#didn't lose any mark status info
        G.unk=NA
        y.sight.unk=NA
        IDunk=NA
      }
    }else if(pMarkID[2]==1&pMarkID[1]<1){#Can't perfectly ID mark status of marked but can unmarked
      stop("It seems unrealistic that you can perfectly determine the marked status unmarked individuals, but not marked individuals. Did not program this possibility.")
    }else{#can perfectly ID mark status of all guys
      y.sight.unk=NA
      G.unk=NA
      IDunk=NA
    }
    #Lose ID's of marked status guys
    if(pID<1){
      idx1=which(y.sight.marked>0)#used to extract counts
      count=y.sight.marked[idx1]
      idx2=which(y.sight.marked>0,arr.ind=TRUE)#used to move counts
      idx3=rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2=idx2[idx3,]
      rem=1-rbinom(nrow(idx2),1,pID[1])#which counts can't we tell marked status
      if(sum(rem)>0){
        idx2=idx2[which(rem==1),]
        if(!is.matrix(idx2)){
          idx2=matrix(idx2,nrow=1)
        }
        IDmnoID=idx2[,1]
        y.sight.marked.noID=array(0,dim=c(nrow(idx2),J,K))
        G.marked.noID=matrix(NA,nrow=nrow(idx2),ncol=ncat)
        for(i in 1:nrow(idx2)){
          #remove no ID guys from known ID marked sightings
          y.sight.marked[idx2[i,1],idx2[i,2],idx2[i,3]]=y.sight.marked[idx2[i,1],idx2[i,2],idx2[i,3]]-1
          #add no ID marked guys to new
          y.sight.marked.noID[i,idx2[i,2],idx2[i,3]]=1
          #add no ID marked guy IDcovs to new structure
          G.marked.noID[i,]=G.marked[idx2[i,1],]
        }
      }else{
        warning("pID<1, but no marked IDs lost during simulation")
        y.sight.marked.noID=NA
        G.marked.noID=NA
        IDmnoID=NA
      }
    }else{
      y.sight.marked.noID=NA
      G.marked.noID=NA
      IDmnoID=NA
    }
    #Amplification failure in unmarked, unknown, and marked no ID if present
    for(l in 1:ncat){
      G.unmarked[which(rbinom(n.samples,1,pIDcat[l])==0),l]=0 #0 is dropout
    }
    if(!is.na(G.unk[1])){
      for(l in 1:ncat){
        G.unk[which(rbinom(nrow(G.unk),1,pIDcat[l])==0),l]=0 #0 is dropout
      }
    }
    if(!is.na(G.marked.noID[1])){
      for(l in 1:ncat){
        G.marked.noID[which(rbinom(nrow(G.marked.noID),1,pIDcat[l])==0),l]=0 #0 is dropout
      }
    }
    #Telemetry observations
    if(tlocs>0){
      locs=array(NA,dim=c(n.marked,tlocs,2))
      for(i in 1:n.marked){
        for(j in 1:tlocs){
          locs[i,j,]=c(rnorm(1,s[IDmarked[i],1],sigma),rnorm(1,s[IDmarked[i],2],sigma))
        }
      }
    }else{
      locs=NA
    }
    
    #check data
    y.sight.check=y.sight*0
    for(i in 1:length(IDmarked)){
      y.sight.check[IDmarked[i],,]=y.sight.marked[i,,]
    }
    if(length(IDum)>0){
      for(i in 1:length(IDum)){
        y.sight.check[IDum[i],,]=y.sight.check[IDum[i],,]+y.sight.unmarked[i,,]
      }
    }
    if(pMarkID[2]<1|pMarkID[1]<1){
      for(i in 1:length(IDunk)){
        y.sight.check[IDunk[i],,]=y.sight.check[IDunk[i],,]+y.sight.unk[i,,]
      }
    }
    if(pID<1){
      for(i in 1:length(IDmnoID)){
        y.sight.check[IDmnoID[i],,]=y.sight.check[IDmnoID[i],,]+y.sight.marked.noID[i,,]
      }
    }
    if(!all(y.sight==y.sight.check)){
      stop("Error rebuilding data. Report bug.")
    }
   
    
    
    out<-list(y.sight=y.sight,y.sight.marked=y.sight.marked,y.sight.unmarked=y.sight.unmarked,
              y.sight.unk=y.sight.unk,y.sight.marked.noID=y.sight.marked.noID,
              G.marked=G.marked,G.unmarked=G.unmarked,
              G.unk=G.unk,G.marked.noID=G.marked.noID,
              IDlist=list(ncat=ncat,IDcovs=IDcovs),locs=locs,
              X=X,K=K,buff=buff,obstype=obstype,s=s,
              IDmarked=IDmarked,IDum=IDum,IDunk=IDunk,IDmnoID=IDmnoID,n.marked=n.marked)
    return(out)
  }