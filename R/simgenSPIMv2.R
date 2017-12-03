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
#' @param theta_c the probability that all of the samples from an individual will be perfectly identified
#' @param theta_p the partial identity thinning rate
#' @param theta_m the proportion of the identity constraint matrix that will have entry 1, allowing subcluster i and j to possibly
#' be the same individual
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from a SCR study with potentially partial individual identities. If theta_c=1, normal SCR data
#' are produced, but if theta_c<1, some individuals will only be partially identified.  If theta_c=0, pgroup=0, and theta_m=0, unidentified
#'  counts or detections (depending on observation model) will be produced as treated by Chandler and Royle (2013).  If pgroup>0, some
#'  singleton samples will be grouped back together.  This reproduces the scenario where you know some samples are the same individual,
#'  but other samples may belong to this focal individual as well.  We can also rule out matches between samples using a constraint
#'  matrix with elements 1 allowing matches and 0 disallowing them.  Because all the samples from the same individual must be allowed to 
#'  match with each other, a certain proportion of 1's are mandatory.  If theta_m is greater than this proportion, random matches
#'  between subclusters are added (1's in the identity constraint matrix) until theta_m percent of the identity constraint matrix
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
  function(N=50,lam0=0.2,sigma=0.50,K=10,X=X,buff=2,theta_c=0.5,theta_p=0.3,maxC=10,theta_m=0.3,obstype="poisson"){
    library(abind)
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    # Simulate encounter history
    y.true <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.true[i,j,k]=rbinom(1,1,pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y.true[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    #remove uncaptured guys
    s.status=rep(NA,N)
    if(any(rowSums(y.true)==0)){
      s.status[rowSums(y.true)==0]="uncaptured"
      idx=which(rowSums(y.true)>0)
      y.true=y.true[idx,,]
    }
    #Figure out who can be identified fully
    if(theta_c>0){
      IDd=rbinom(nrow(y.true),1,theta_c)
      y.complete=y.true[IDd==1,,]
      y.partial=y.true[IDd==0,,]
      s.status[idx[IDd==1]]="complete"
      s.status[idx[IDd==0]]="partial"
      if(length(dim(y.partial))==2){
        y.partial=array(y.partial,dim=c(1,J,K))
      }
    }else{
      y.compelete=array(0,dim=c(0,K,J))
      y.partial=y.true
    }
    if(any(rowSums(y.partial)>0)){
      y.partial=y.partial[rowSums(y.partial)>0,,]
      if(length(dim(y.partial))==2){
        y.partial=array(y.partial,dim=c(1,J,K))
      }
    }
    #If any partial samples
    if(sum(y.partial)>0){
      #Break identity linkeages
      if(theta_p>0){
        #zero-truncated binomial simulator
        rtbin <- function(n, N, p) {
          if(N==0)
            return(1)
          X0 <- 0:N
          TB0 <- dbinom(X0, N, p)[-1]
          TB <- TB0/sum(TB0)
          sample(1:N, size=n, replace=TRUE, prob=TB)
        }
        y.partial.true <- array(0L, dim=c(N, J, K, maxC))
        y.partial.data.list <- list()
        N.partial=nrow(y.partial)
        detected <- logical(N.partial)
        nDets <- integer(N.partial)
        C <- rep(1L, N.partial)
        ## Binary cluster indicator
        c <- matrix(0L, N.partial, maxC)
        det.ctr <- 1
        id <- integer(0)
        for(i in 1:N.partial) {
          y.partial.i <- y.partial[i,,]
          C[i] <- rtbin(1, sum(y.partial.i), theta_p) ## nClusters
          if(C[i] >= maxC){
            warning("Need to raise maxC")
          }
          c[i,1:C[i]] <- 1
          for(j in 1:J) {
            for(k in 1:K) {
              if(y.partial.i[j,k]<1)
                next
              y.partial.true[i,j,k,] <- rmultinom(1, y.partial.i[j,k], c[i,])
            }
          }
          for(l in 1:maxC) {
            if(all(y.partial.true[i,,,l]<1))
              next
            y.partial.data.list[[det.ctr]] <- y.partial.true[i,,,l]
            id[det.ctr] <- i
            det.ctr <- det.ctr+1
          }
        }
        clusterCounts <- apply(y.partial.true, c(1,4), sum)
        nDets <- rowSums(y.partial)
        L <- length(y.partial.data.list)
        y.partial.data <- array(NA, c(L, J, K))
        for(l in 1:L) {
          y.partial.data[l,,] <- y.partial.data.list[[l]]
        }
      }else{#No fissions of partial ID guys
        y.partial.data=y.partial
      }
      #Add identity constraints
      N.partial.data=nrow(y.partial.data)
      #default is only partials from same ID can match
      constraints=matrix(1,nrow=N.partial.data,ncol=N.partial.data)
      for(i in 1:N.partial.data){
        constraints[i,id!=id[i]]=0
      }
      #Add possiblity of matching between different IDs
      if(theta_m>0){
        nconstrain=round(theta_m*N.partial.data^2)
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
          warning("theta_m too small to produce number of constraints more than complete identities imply")
        }
      }
    }else{
      y.partial.data=y.partial
      id=NA
      constraints=NA
    }
    #Fix dimension of y.complete if no complete guys
    if(length(dim(y.complete))==2){
      y.complete=array(y.complete,dim=c(1,dim(y.complete)))
    }
    out<-list(y.true=y.true,y.complete=y.complete,y.partial=y.partial,y.partial.data=y.partial.data,constraints=constraints,s=s,s.status=s.status,X=X, K=K,buff=buff,
              obstype=obstype,ID=id)
    return(out)
  }

