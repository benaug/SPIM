e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#' Simulate data from a SCR study
#' @param N a vector indicating the number of individuals to simulate
#' @param lam0 the detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the K x 2 matrix of trap locations
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from a camera trap SCR study. The extent of the state space is controlled by "buff", which buffers the
#' minimum and maximum X and Y extents.  Therefore, it probably only makes sense for square or rectangular grids.  Functionality
#' for user-provided polygonal state spaces will be added in the future.
#' @author Ben Augustine
#' @export

simgenSPIM <-
  function(N=50,lam0=0.2,sigma=0.50,K=10,X=X,buff=2,pID=0.5,pconstrain=0.3,obstype="bernoulli"){
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
    if(pID>0){
      IDd=rbinom(N,1,pID)
      yID=y[IDd==1,,]
      yUnk=y[IDd==0,,]
      if(length(dim(yUnk))==2){
        yUnk=array(yUnk,dim=c(1,J,K))
      }
    }else{
      yID=array(0,dim=c(0,K,J))
      yUnk=y
    }
    if(any(rowSums(yID)>0)){
      yID=yID[rowSums(yID)>0,,]
    }
    if(any(rowSums(yUnk)>0)){
      yUnk=yUnk[rowSums(yUnk)>0,,]
      if(length(dim(yUnk))==2){
        yUnk=array(yUnk,dim=c(1,J,K))
      }
    }
    if(sum(yUnk)>0){
      caps=which(yUnk==1,arr.ind=TRUE)
      nobs=nrow(caps)
      yUnkobs=array(0,dim=c(nobs,J,K))
      for(i in 1:nobs){
        yUnkobs[i,caps[i,2],caps[i,3]]=1
      }
      #Find constraints based on all known identities
      constraints=matrix(1,nrow=nobs,ncol=nobs)
      for(i in 1:nobs){
        constraints[i,which(caps[,1]!=caps[i,1])]=0
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
    }
    out<-list(yID=yID,yUnk=yUnk,yUnkobs=yUnkobs,constraints=constraints,s=s,X=X, K=K,buff=buff,obstype=obstype)
    return(out)
  }
