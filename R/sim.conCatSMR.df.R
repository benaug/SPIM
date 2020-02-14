#' Simulate data from the conventional (no marking process) categorical spatial mark resight model with or 
#' without telemetry data and with a detection function that varies by one categorical
#' identity covariate. See sim.conCatSMR() for the same function with a single detection function.
#' @param N Abundance
#' @param n.marked The number of marked individuals in the population, distributed randomly across the state space.
#' an error will be produced if fewer than n.marked individuals are captured and marktype=="natural".
#' @param lam0 a vector of detection function baseline detection rates the same length as the first list 
#' element in IDcovs. Converted to p0, the baseline detection probability if the obstype="bernoulli".
#' @param sigma a vector of detection function spatial scale parameters the same length as the first list 
#' element in IDcovs.
#' @param K The number of sampling occasions
#' @param X a matrix with two columns for the X and Y trap locations. J rows.
#' @param buff an integer indicating the distance to buffer the trapping array, X, to create the state space
#' @param obstype a character string indicating the observation model "bernoulli" or "poisson"
#' @param ncat an integer indicating the number of categorical identity covariates
#' @param pIDcat a vector of length ncat containing the probability that the value of each categorical identity covariate
#'  is observed upon capture
#' @param IDcovs a list of length ncat containing the values each categorical identity covariate can take. The
#' length of each list element determines the number of values each categorical identity covariate can take.
#' @param gamma a list of the category level probabilities of the same dimension as IDcovs. The category level
#' probabilities for each covariate must sum to 1.
#' @param pMarkID a vector of length 2 containing the probability the marked status is observed upon capture
#' for marked and unmarked individuals, respectively. If these are less than 1, unknown marked status
#' samples are produced.
#' @param tlocs a single integer indicating the number of telemetry locations to simulate for each marked individual.
#' @param pID the probability a marked individual's identity is obtained upon capture. If this is less
#' than one, marked but unknown identity samples are produced.
#' @param marktype a character string indicating whether marks are preallocated or obtained from natural marks, "premarked", or "natural".
#' @description This function simulates data from a conventional spatial mark resight survey for categorically
#' marked populations where the detection function parameters vary by the levels of the first categorical
#' identity covariate. Whether marks are preallocated or identified upon capture (natural
#' marks) is determined through "marktype". Imperfect determination of marked status is controlled through "pMarkID"
#' and imperfect individual identification of marked individuals is controlled through "pID". Telemetry data for
#' marked indiviuals is added through "tlocs".
#' @return a list with many elements. y.sight is the complete sighting history for all individuals. 
#' y.sight.marked is the sighting history of the marked, observed marked status, and individually identified samples.
#' y.sight.unmarked is the sighting history of the observed marked status unmarked individual samples.
#' y.sight.unk is the sighting history for the samples for which marked status could not be determined.
#' y.sight.marked.noID is the sighting history of the observed marked status marked, but not individually identified samples.
#' Not all structures will be produced if there is perfect observation of mark status and/or individual identity of 
#' marked individuals. y.sight is of dimension N x J x K. y.sight.marked si of dimension
#' n.marked X J x K. y.sight.unmarked is of dimension n_um x J x K, 
#' where n_um is the number of unmarked samples observed. The other latent identity sighting histories
#' are similarly structured with one observation per i.
#' 
#' G.x structures housing the observed categorical identity covariates correspond to the y.sight.x 
#' structures, linked by the i dimension. "s" contains the simulated activity centers corresponding 
#' to y.sight, the complete, perfect identity data. Missing values, if simulated, are indicated with a 0.
#' 
#' IDlist is a list containing ncat and IDcovs, inputs to the simulation function.
#' 
#' IDmarked, IDum, IDunk, and IDmnoID indicate which individual
#' in "s" each ith row of the latent identity sighting histories came from. These could be used to
#' reassemble the latent identity data sets into y.sight.
#' 
#' "locs" contains an n.marked x nlocs x 2 array of telemetry locations, if simulated, for the marked
#' individuals. The i dimension of locs corresponds to the first n.marked i elements of y.sight and the
#' i dimension of y.sight.marked.
#' @author Ben Augustine
#' @examples
#' \dontrun{
#' N=100
#' n.marked=40
#' lam0=c(0.1,0.3) #Enter same number of detection function parameters as nlevels[i]
#' sigma=c(0.7,0.5) #same number of lam0 and sigma parameters
#' K=20 #occasions
#' buff=3 #state space buffer
#' X<- expand.grid(3:11,3:11) #trapping array
#' pMarkID=c(1,1) #probability of observing marked status of marked and unmarked individuals
#' pID=1 #Probability marked individuals are identified
#' ncat=3  #number of ID categories
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each category level. Assume equal here.
#' nlevels=rep(2,ncat) #number of levels per IDcovs
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#category observation probabilities
#' tlocs=0 #telemetry locs/marked individual
#' obstype="poisson" #observation model, count or presence/absence?
#' marktype="premarked" #premarked or natural ID (marked individuals must be captured)?
#' data=sim.conCatSMR.df(N=N,n.marked=n.marked,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,ncat=ncat,
#'                          pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs,pID=pID,marktype=marktype)
#' }
#' @export
sim.conCatSMR.df <-
  function(N=50,n.marked=10,lam0=0.25,sigma=0.50,K=10,X=X,buff=3,obstype="bernoulli",ncat=ncat,
           pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=c(1,1),tlocs=0,pID=1,
           marktype="premarked"){
    if(!marktype%in%c("natural","premarked")){
      stop("marktype must be 'natural' or 'premarked'")
    }
    #simulate IDcat
    G.true=matrix(NA,nrow=N,ncol=ncat) #all IDcat in population.
    for(i in 1:N){
      for(j in 1:ncat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }
    library(abind)
    # simulate a population of activity centers
    X=as.matrix(X)
    J=nrow(X)
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    #Variable df stuff
    if(length(lam0)!=length(sigma)){
      stop("Need the same number of values for each detection function parameter")
    }
    if(length(IDcovs[[1]])!=length(lam0)){
      stop("IDcovs[[1]] should have the same number of values as detection function parameter values")
    }
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    for(i in 1:length(sigma)){
      lamd[G.true[,1]==i,]<- lam0[i]*exp(-D[G.true[,1]==i,]^2/(2*sigma[i]*sigma[i]))
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
      if(length(idx)<n.marked){
        stop("Fewer than n.marked individuals captured. Cannot naturally mark uncaptured individuals.")
      }
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
                G.unmarked[idx,]=G.true[umguys[i],]
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
        if(length(dim(y.sight.unk))==2){
          y.sight.unk=array(y.sight.unk,dim=c(1,J,K))
        }
        if(length(dim(y.sight.unmarked))==2){
          y.sight.unmarked=array(y.sight.unmarked,dim=c(1,J,K))
        }
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
        if(ncol(G.unmarked)!=ncat){
          G.unmarked=t(G.unmarked)
        }
        if(ncol(G.unk)!=ncat){
          G.unk=t(G.unk)
        }
        y.sight.unk=y.sight.unmarked[unm2unk,,]
        y.sight.unmarked=y.sight.unmarked[-unm2unk,,]
        if(length(dim(y.sight.unk))==2){
          y.sight.unk=array(y.sight.unk,dim=c(1,J,K))
        }
        if(length(dim(y.sight.unmarked))==2){
          y.sight.unmarked=array(y.sight.unmarked,dim=c(1,J,K))
        }
        IDunk=IDum[unm2unk]
        IDum=IDum[-unm2unk]
      }else{#didn't lose any mark status info
        y.sight.unk=NA
        G.unk=NA
        IDunk=NA
      }
      #add marked guys
      if(sum(y.sight.marked)==0)next
      idx1=which(y.sight.marked>0)#used to extract counts
      count=y.sight.marked[idx1]
      idx2=which(y.sight.marked>0,arr.ind=TRUE)#used to move counts
      idx3=rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2=idx2[idx3,]
      if(!is.matrix(idx2)){
        idx2=matrix(idx2,ncol=3)
      }
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
        if(ncol(G.unk2)!=ncat){
          G.unk2=t(G.unk2)
        }
        if(all(!is.na(y.sight.unk))&all(!is.na(y.sight.unk2))){
          y.sight.unk=abind(y.sight.unk,y.sight.unk2,along=1)
          G.unk=rbind(G.unk,G.unk2)
          IDunk=c(IDunk,idx2[,1])
        }else if(!is.null(y.sight.unk2)){
          y.sight.unk=y.sight.unk2
          G.unk=G.unk2
          IDunk=idx2[,1]
        }
      }
    }else if(pMarkID[2]==1&pMarkID[1]<1){#Can't perfectly ID mark status of marked but can unmarked
      stop("It seems unrealistic that you can perfectly determine the marked status unmarked individuals, but not marked individuals. Did not program this possibility.")
    }else{#can perfectly ID mark status of all guys
      y.sight.unk=NA
      G.unk=NA
      IDunk=NA
    }
    #Lose ID's of marked status guys
    if(pID<1&sum(y.sight.marked>0)){
      idx1=which(y.sight.marked>0)#used to extract counts
      count=y.sight.marked[idx1]
      idx2=which(y.sight.marked>0,arr.ind=TRUE)#used to move counts
      idx3=rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2=idx2[idx3,]
      if(!is.matrix(idx2)){
        idx2=matrix(idx2,ncol=3)
      }
      if(nrow(idx2)>0){
        rem=1-rbinom(nrow(idx2),1,pID[1])#which counts can't we tell marked status
      }else{
        rem=0
      }
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
      drop=which(rbinom(nrow(G.unmarked),1,pIDcat[l])==0)
      if(length(drop)>0){
        G.unmarked[drop,l]=0 #0 is dropout
      }
    }
    if(!is.na(G.unk[1])){
      for(l in 1:ncat){
        drop=which(rbinom(nrow(G.unk),1,pIDcat[l])==0)
        if(length(drop)>0){
          G.unk[drop,l]=0 #0 is dropout
        }
      }
    }
    if(!is.na(G.marked.noID[1])){
      for(l in 1:ncat){
        drop=which(rbinom(nrow(G.marked.noID),1,pIDcat[l])==0)
        if(length(drop)>0){
          G.marked.noID[drop,l]=0 #0 is dropout
        }
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
    if(all(!is.na(IDunk))){
      for(i in 1:length(IDunk)){
        y.sight.check[IDunk[i],,]=y.sight.check[IDunk[i],,]+y.sight.unk[i,,]
      }
    }
    if(pID<1){
      if(all(is.finite(IDmnoID))){
        for(i in 1:length(IDmnoID)){
          y.sight.check[IDmnoID[i],,]=y.sight.check[IDmnoID[i],,]+y.sight.marked.noID[i,,]
        }
      }
    }
    if(!all(y.sight==y.sight.check)){
      stop("Error rebuilding data. Report bug.")
    }
   
    
    
    out<-list(y.sight=y.sight,y.sight.marked=y.sight.marked,y.sight.unmarked=y.sight.unmarked,
              y.sight.unk=y.sight.unk,y.sight.marked.noID=y.sight.marked.noID,G.true=G.true,
              G.marked=G.marked,G.unmarked=G.unmarked,
              G.unk=G.unk,G.marked.noID=G.marked.noID,
              IDlist=list(ncat=ncat,IDcovs=IDcovs),locs=locs,
              X=X,K=K,buff=buff,obstype=obstype,s=s,
              IDmarked=IDmarked,IDum=IDum,IDunk=IDunk,IDmnoID=IDmnoID,n.marked=n.marked)
    return(out)
  }