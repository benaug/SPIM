e2dist<-function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#' Simulate data with potentially partial individual indentities
#' @param N a vector indicating the number of individuals to simulate
#' @param lam0 the detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the K x 2 matrix of trap locations
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @param pID the probability that all of the samples from an individual will be perfectly identified
#' @param pgroup the proportion of partial identity samples that will be randomly grouped together in to correct identity subclusters
#' @param pconstrain the proportion of the identity constraint matrix that will have entry 1, allowing subcluster i and j to possibly
#' be the same individual
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from a SCR study with potentially partial individual identities. If pID=1, normal SCR data
#' are produced, but if pID<1, some individuals will only be partially identified.  If pID=0, pgroup=0, and pconstrain=0, unidentified
#'  counts or detections (depending on observation model) will be produced as treated by Chandler and Royle (2013).  If pgroup>0, some
#'  singleton samples will be grouped back together.  This reproduces the scenario where you know some samples are the same individual,
#'  but other samples may belong to this focal individual as well.  We can also rule out matches between samples using a constraint
#'  matrix with elements 1 allowing matches and 0 disallowing them.  Because all the samples from the same individual must be allowed to 
#'  match with each other, a certain proportion of 1's are mandatory.  If pconstrain is greater than this proportion, random matches
#'  between subclusters are added (1's in the identity constraint matrix) until pconstrain percent of the identity constraint matrix
#'  are 1s.
#'  
#' The extent of the state space is controlled by "buff", which buffers the minimum and maximum X and Y extents.  Therefore, 
#' it probably only makes sense for square or rectangular grids.  Functionality for user-provided polygonal state spaces will 
#' be added in the future.
#' 
#' While the binomial observation model is included, I'm not sure it makes sense in practice unless it can be ensured
#' that multiple capture events from the same individual at the same trap on the same occasion never happen.  If you 
#' cannot fully identify all samples, you can't convert all trap-occasion counts>1 to a single encounter.
#' 
#' @author Ben Augustine
#' @export

simgenSPIM <-
  function(N=50,lam0=0.2,sigma=0.50,K=10,X=X,buff=2,pID=0.5,theta=0.3,pconstrain=0.3,obstype="poisson"){
    library(abind)
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    # Simulate encounter history
    y <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,j,k]=rbinom(1,1,pd[i,j])
            }
          }
        }
    }else if(obstype=="poisson"){
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,j,k]=rpois(1,lamd[i,j])
            }
          }
        }
    }else{
      stop("observation model not recognized")
    }
    #Figure out who can be identified fully
    s.status=rep(NA,N)
    if(any(rowSums(y)==0)){
      s.status[rowSums(y)==0]="uncaptured"
      idx=which(rowSums(y)>0)
      y=y[idx,,]
    }
    if(pID>0){
      IDd=rbinom(nrow(y),1,pID)
      yID=y[IDd==1,,]
      yUnk=y[IDd==0,,]
      s.status[idx[IDd==1]]="ID"
      s.status[idx[IDd==0]]="unk"
      if(length(dim(yUnk))==2){
        yUnk=array(yUnk,dim=c(1,J,K))
      }
    }else{
      yID=array(0,dim=c(0,K,J))
      yUnk=y
    }
    if(any(rowSums(yUnk)>0)){
      yUnk=yUnk[rowSums(yUnk)>0,,]
      if(length(dim(yUnk))==2){
        yUnk=array(yUnk,dim=c(1,J,K))
      }
    }
    if(sum(yUnk)>0){
      caps=which(yUnk>0,arr.ind=TRUE)
      unkIDs=c()
      yUnkobs=array(0,dim=c(0,J,K))
      if(theta>0){#Break apart yUnk
        idx2=1
        for(i in 1:nrow(yUnk)){
          maxbreak=sum(yUnk[i,,])-1#number of possible breaks
          B=rbinom(1,maxbreak,theta)#number of breaks
          if(B>0){
            # deal=maxbreak+1-(B+1)
            # C=rmultinom(1,deal,rep(1/(B+1),B+1))
            guys=which(yUnk[i,,]>0)#indices with counts
            guys=rep(guys,times=yUnk[i,,][guys])#break apart poisson counts, too
            if(length(guys)>1){
              guys2=sample(guys,length(guys))#shuffled sample numbers
            }else{
              guys2=guys
            }
            # idx=rep(1:(B+1),times=C+1)
          
            #make an index to break guys into B+1 groups
            if(length(guys2)==(B+1)){#if there are as many samples as groups, sample without replacement
              idx=sample(1:(B+1),length(guys2),replace=FALSE)
            }else{#otherwise, need to make sure all samples occur at least once, then sample with replacement
              idx=c(1:(B+1),sample(1:(B+1),length(guys2)-B-1,replace=TRUE))#make sure all members of 1:(B+1) occur
              idx=idx[sample(1:length(idx),length(idx))]#shuffle this vector before selecting index
            }
            
            for(j in 1:(B+1)){#use idx to subset out into B+1 groups
              yUnkobs=abind(yUnkobs,matrix(0,nrow=J,ncol=K),along=1)
              # yUnkobs[idx2,,][guys2[idx==j]]=1
              guys3=guys2[idx==j]
              for(k in 1:length(idx)){
                yUnkobs[idx2,,][guys3[k]]=yUnkobs[idx2,,][guys3[k]]+1

              }
              idx2=idx2+1
              unkIDs=c(unkIDs,i)
            }
          }else{
            yUnkobs=abind(yUnkobs,matrix(0,nrow=J,ncol=K),along=1)
            yUnkobs[idx2,,]=yUnk[i,,]
            idx2=idx2+1
            unkIDs=c(unkIDs,i)
          }
        }
      }else{
        unkIDs=1:nrow(yUnk)
        yUnkobs=yUnk
      }
      
     
      #Find constraints based on all known identities
      nobs=nrow(yUnkobs)#changes if some samples associated
      constraints=matrix(1,nrow=nobs,ncol=nobs)
      for(i in 1:nobs){
        constraints[i,unkIDs!=unkIDs[i]]=0
      }
      if(pconstrain>0){
        nconstrain=round(pconstrain*nobs^2)
        obsconstrain=sum(constraints)
        #Add some extra possible matches between different individuals
        if(nconstrain>obsconstrain){
          deal=round((nconstrain-obsconstrain)/2) #deal half to upper triangle matrix and fill in respective lower traingle matrix
          zeros=which(constraints==0&upper.tri(constraints),arr.ind=TRUE)
          choose=sample(1:nrow(zeros),deal)
          for(i in 1:deal){
            constraints[zeros[choose[i],1],zeros[choose[i],2]]=1 #upper tri
            constraints[zeros[choose[i],2],zeros[choose[i],1]]=1 #lower tri
          }
        }else{
          warning("pconstrain too small to produce number of constraints more than complete identities imply")
        }
      }
    }else{
      warning("No unknown identity individuals produced")
      constraints=NA
      yUnkobs=yUnk
      unkIDs=NA
    }
    if(length(dim(yID))==2){
      yID=array(yID,dim=c(1,dim(yID)))
    }
    out<-list(yID=yID,yUnk=yUnk,yUnkobs=yUnkobs,constraints=constraints,s=s,s.status=s.status,X=X, K=K,buff=buff,
              obstype=obstype,unkIDs=unkIDs)
    return(out)
  }

