#' Simulate data from the generalized categorical spatial mark resight model with or 
#' without telemetry data
#' @param N abundance
#' @param lam0.mark the detection function baseline detection rate for the marking process. Converted to p0, the baseline detection probability if the obstype="bernoulli".
#' @param lam0.sight the detection function baseline detection rate for the sighting process. Converted to p0, the baseline detection probability if the obstype="bernoulli".
#' @param sigma the detection function spatial scale parameter assumed to be the same for both capture
#' processes. This could be relaxed.
#' @param K1 the number of marking occasions
#' @param K2 the number of sighting occasions
#' @param Korder a character vector of length K1+K2 specifying the order of the marking and sighting occasions.
#' Vector elements are either "M" or "S" arranged in the order the marking and sighting occasions occurred.
#' See example below. Korder does not need to be an input if all sighting occasions occurred after the final 
#' marking occasion.
#' @param X1 a matrix with two columns for the X and Y trap locations of the marking process. J1 rows.
#' @param X2 a matrix with two columns for the X and Y trap locations of the sighting process. J2 rows.
#' @param buff an integer indicating the distance to buffer the combined trapping array (X1 and X2), to create the state space
#' @param obstype a vector of length two indicating the observation model, "bernoulli" or "poisson", for the 
#' marking and sighting process
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
#' @description This function simulates data from a generalized spatial mark resight survey for categorically
#' marked populations. Capture histories from the marking and sighting process are produced. If there is
#' only one categorical identity covariate with one value, the function simulates
#' from typical generalized mark resight. Imperfect determination of marked status is controlled through "pMarkID"
#' and imperfect individual identification of marked individuals is controlled through "pID". Telemetry data for
#' marked indiviuals is added through "tlocs".
#' @return a list with many elements. y.mark is the capture history for the marking process where all
#' individual identities of captured individuals are known or unidentified individuals are excluded.
#' y.mark is the capture history from the marking process.  
#' y.sight is the complete sighting history for all captured individuals. 
#' y.sight.marked is the sighting history of the marked, observed marked status, and individually identified samples.
#' y.sight.unmarked is the sighting history of the observed marked status unmarked individual samples.
#' y.sight.unk is the sighting history for the samples for which marked status could not be determined.
#' y.sight.marked.noID is the sighting history of the observed marked status marked, but not individually identified samples.
#' Not all structures will be produced if there is perfect observation of mark status and/or individual identity of 
#' marked individuals. y.mark is of dimension n.marked x J1 x K1, where n.marked is the number of 
#' individuals captured and marked. y.sight is of dimension n.total x J2 x K2, where n.total is the total
#' number of individuals captured and sighted. y.sight.unmarked is of dimension n_um x J2 x K2, 
#' where n_um is the number of unmarked samples observed. The other latent identity sighting histories
#' are similarly structured with one observation per i.
#' 
#' G.x structures housing the observed categorical identity covariates correspond to the y.sight.x 
#' structures, linked by the i dimension. Missing values, if simulated, are indicated with a 0.
#' 
#' IDlist is a list containing ncat and IDcovs, inputs to the simulation function.
#' 
#' IDmarked, IDum, IDunk, and IDmnoID indicate which individual
#' in "sfull" each ith row of the latent identity sighting histories came from. These could be used to
#' reassemble the latent identity data sets into y.sight.
#' 
#' "locs" contains an n.marked x nlocs x 2 array of telemetry locations, if simulated, for the marked
#' individuals. The i dimension of locs corresponds to the first n.marked i elements of y.sight and the
#' i dimension of y.sight.marked.
#' 
#' "markedS" contains a history of which of the marked individuals were marked on each occasion.
#' @author Ben Augustine
#' @examples
#' \dontrun{
#' N=50
#' lam0.mark=0.05 #marking process baseline detection rate
#' lam0.sight=0.25 #sighting process baseline detection rate
#' sigma=0.50 #shared spatial scale parameter
#' K1=10 #number of marking occasions
#' K2=10 #number of sighting occasions
#' buff=2 #state space buffer
#' X1<- expand.grid(3:11,3:11) #marking locations
#' X2<- expand.grid(3:11+0.5,3:11+0.5) #sighting locations
#' pMarkID=c(0.8,0.8) #probability of observing marked status of marked and unmarked individuals
#' pID=0.8 #probability of determining identity of marked individuals
#' obstype=c("bernoulli","poisson") #observation model of both processes
#' ncat=3  #number of categorical identity covariates
#' gamma=IDcovs=vector("list",ncat) 
#' nlevels=rep(2,ncat) #number of IDcovs per loci
#' for(i in 1:ncat){ 
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' #inspect ID covariates and level probabilities
#' str(IDcovs) #3 covariates with 2 levels
#' str(gamma) #each of the two levels are equally probable
#' pIDcat=rep(1,ncat) #category observation probabilities
#' #Example of interspersed marking and sighting. 
#' Korder=c("M","M","S","S","S","S","M","M","S","M","M","S","M","M","S","S","S","S","M","M")
#' #Example with no interspersed marking and sighting.
#' Korder=c(rep("M",10),rep("S",10))
#' #if Korder is not of length K1+K2, there will be an error.
#' #Omitting Korder assumes no sighting occurred before marking ended.
#' tlocs=25
#' data=sim.genCatSMR(N=N,lam0.mark=lam0.mark,lam0.sight=lam0.sight,sigma=sigma,K1=K1,
#'                    K2=K2,Korder=Korder,X1=X1,X2=X2,buff=buff,obstype=obstype,ncat=ncat,
#'                    pIDcat=pIDcat,IDcovs=IDcovs,gamma=gamma,pMarkID=pMarkID,pID=pID,tlocs=tlocs)
#'data$markedS #marked status history across occasions for marked individuals
#'}
#' @export

sim.genCatSMR <-
  function(N=50,lam0.mark=0.25,lam0.sight=0.25,sigma=0.50,K1=10,K2=10,Korder=NA,X1=X1,X2=X2,buff=3,obstype="bernoulli",ncat=ncat,
           pIDcat=pIDcat,IDcovs=IDcovs,gamma=gamma,pMarkID=c(1,1),tlocs=0,pID=1){
    library(abind)
    # simulate a population of activity centers
    X1=as.matrix(X1)
    X2=as.matrix(X2)
    Xall=rbind(X1,X2)
    s<- cbind(runif(N, min(Xall[,1])-buff,max(Xall[,1])+buff), runif(N,min(Xall[,2])-buff,max(Xall[,2])+buff))
    D1<- e2dist(s,X1)
    D2<- e2dist(s,X2)
    J1=nrow(X1)
    J2=nrow(X2)
    if(is.na(Korder[1])){
      warning("Assuming all marking occasions preceded the first sighting occasion since 'Korder' was omitted")
      Korder=c(rep("M",K1),rep("S",K2))
    }
    if(!all(c("M","S")%in%names(table(Korder)))|!all(names(table(Korder))%in%c("M","S"))){
      stop("Korder must only contain characters'M'and 'S' indicating marking and sighting sessions")
    }
    if(length(Korder)!=sum(K1,K2)){
      stop("Korder is not the right size")
    }
    if((sum(Korder=="M")!=K1)|(sum(Korder=="S")!=K2)){
      stop("Fix number of M and S in Korder to match K1 and K2")
    }
    #simulate IDcovs
    G.true=matrix(NA,nrow=N,ncol=ncat) #all IDcovs in population.
    for(i in 1:N){
      for(j in 1:ncat){
        G.true[i,j]=sample(IDcovs[[j]],1,prob=gamma[[j]])
      }
    }

    for(i in 1:length(sigma)){
      lamd.mark<- lam0.mark*exp(-D1^2/(2*sigma*sigma))
      lamd.sight<- lam0.sight*exp(-D2^2/(2*sigma*sigma))
    }
    # Capture and mark individuals
    y.mark <-array(0,dim=c(N,J1,K1))
    if(obstype[1]=="bernoulli"){
      pd.mark=1-exp(-lamd.mark)
      for(i in 1:N){
        for(j in 1:J1){
          for(k in 1:K1){
            y.mark[i,j,k]=rbinom(1,1,pd.mark[i,j]) 
          }
        }
      }
    }else if(obstype[1]=="poisson"){
      for(i in 1:N){
        for(j in 1:J1){
          for(k in 1:K1){
            y.mark[i,j,k]=rpois(1,lamd.mark[i,j])
          }
        }
      }
    }else{
      stop("obstype 1 not recognized")
    }
    #Sight individuals
    y.sight <-array(0,dim=c(N,J2,K2))
    if(obstype[2]=="bernoulli"){
      pd.sight=1-exp(-lamd.sight)
      for(i in 1:N){
        for(j in 1:J2){
          for(k in 1:K2){
            y.sight[i,j,k]=rbinom(1,1,pd.sight[i,j]) 
          }
        }
      }
    }else if(obstype[2]=="poisson"){
      for(i in 1:N){
        for(j in 1:J2){
          for(k in 1:K2){
            y.sight[i,j,k]=rpois(1,lamd.sight[i,j])
          }
        }
      } 
    }else{
      stop("obstype 2 not recognized")
    }
    #Figure out who is marked on each sighting occasion
    markedT=1*(apply(y.mark,c(1,3),sum)>0)
    for(k in 2:K1){
      markedT[markedT[,k-1]==1,k]=1
    }
    markedS=matrix(0,ncol=K2,nrow=N)
    markocc=which(Korder=="M")
    sightocc=which(Korder=="S")
    for(k in 1:K2){
      prevmark=which(markocc<sightocc[k])
      if(length(prevmark)>0){
        markedS[,k]=markedT[,max(prevmark)]
      }
    }
    #Discard uncaptured individuals (by both methods)
    rem=which(!(rowSums(y.mark)>0|rowSums(y.sight)>0))
    sfull=s
    if(length(rem)>0){
      y.mark=y.mark[-rem,,]
      y.sight=y.sight[-rem,,]
      G.cap=G.true[-rem,]
      markedS=markedS[-rem,]
      s=s[-rem,]
    }
    IDmarked=which(rowSums(y.mark)>0)
    
    #split sightings into marked and unmarked histories, considering occasion of marking
    y.sightMK=apply(y.sight,c(1,3),sum)
    n.marked=sum(rowSums(markedS)!=0)
    n.samples=sum(y.sightMK[(y.sightMK>0&markedS==0)])
    y.sight.marked=0*y.sight
    y.sight.unmarked=array(0,dim=c(n.samples,J2,K2))
    n=nrow(y.sight)
    IDum=rep(NA,n.samples)
    G.unmarked=matrix(NA,nrow=n.samples,ncol=ncat)
    G.marked=matrix(NA,nrow=nrow(y.sight),ncol=ncat)
    G.marked[IDmarked,]=G.true[IDmarked,]
    idx=1
    if(!is.matrix(G.cap)){
      G.cap=matrix(G.cap,ncol=1)
    }
    for(i in 1:nrow(y.sight)){ #loop through inds (uncaptured already removed)
      for(j in 1:J2){ #then traps
        for(k in 1:K2){ #then occasions
          if(y.sight[i,j,k]>0){ #is there at least one sample here?
            if(markedS[i,k]==1){
              G.marked[i,]=G.cap[i,]
              # G.marked[i,]=G.true[i,]
              y.sight.marked[i,j,k]=y.sight[i,j,k]
            }else{
              for(l in 1:y.sight[i,j,k]){ #then samples
                # G.unmarked[idx,]=G.true[i,]
                G.unmarked[idx,]=G.cap[i,]
                y.sight.unmarked[idx,j,k]=1
                IDum[idx]=i
                idx=idx+1
              }
            }
          }
        }
      }
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
    #remove zero histories
    keep=which(rowSums(y.mark)>0)
    y.mark=y.mark[keep,,]
    y.sight.marked=y.sight.marked[keep,,]
    G.marked=G.marked[keep,]
    markedS=markedS[keep,]
    # keep=which(rowSums(y.sight.unmarked)>0)
    # y.sight.unmarked=y.sight.unmarked[keep,,]
    
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
   
    
    
    out<-list(y.mark=y.mark,y.sight=y.sight,y.sight.marked=y.sight.marked,y.sight.unmarked=y.sight.unmarked,
              y.sight.unk=y.sight.unk,y.sight.marked.noID=y.sight.marked.noID,
              G.marked=G.marked,G.unmarked=G.unmarked,
              G.unk=G.unk,G.marked.noID=G.marked.noID,
              IDlist=list(ncat=ncat,IDcovs=IDcovs),locs=locs,
              X1=X1,X2=X2,K1=K1,K2=K2,buff=buff,obstype=obstype,s=s,sfull=sfull,
              IDmarked=IDmarked,
              IDum=IDum,IDunk=IDunk,IDmnoID=IDmnoID,markedS=markedS)
    return(out)
  }