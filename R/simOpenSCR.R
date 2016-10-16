e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#' Simulate data from a Open population SCR study.  Update help file later.  Add polygon SS later.
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

simOpenSCR <-
  function(N=c(40,60,80),gamma=NULL,phi=rep(0.8,2),lam0=rep(0.2,3),sigma=rep(0.50,3),K=rep(10,3),X=X,t=3,M=M,sigma_t=NULL,buff=3,obstype="bernoulli"){
    #Check for user errors
    if(length(N)==1&is.null(gamma)){
      stop("Must provide gamma if length(N)==1")
    }
    if(length(N)==t&!is.null(gamma)){
      stop("Do not provide gamma if length(N)=t")
    }
    storeparms=list(N=N,gamma=gamma,lam0=lam0,sigma=sigma,phi=phi)

    J=unlist(lapply(X,nrow))
    maxJ=max(J)
    #Get state space extent such that buffer is at least buff on all grids
    #minmax=rbind(apply(X[[1]],2,min),apply(X[[1]],2,max))
    minmax=array(NA,dim=c(length(X),2,2))
    for(i in 1:length(X)){
      minmax[i,,]=rbind(apply(X[[i]],2,min),apply(X[[i]],2,max))
    }
    xylim=rbind(apply(minmax,3,min),apply(minmax,3,max))
    xlim=xylim[,1]
    ylim=xylim[,2]
    # # simulate a population of meta ACS and ACs, calculate D and lamd
    mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
    s=array(NA,dim=c(M,t,2))
    D=lamd=array(NA,dim=c(M,maxJ,t))
    if(length(lam0)==1){
      lam0=rep(lam0,t)
    }
    if(length(sigma)==1){
      sigma=rep(sigma,t)
    }
    if(is.null(sigma_t)){#no movement
      for(i in 1:t){
        s[,i,]=mu
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }
    }else{
      for(i in 1:t){#meta mu movement
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,mu[j,1],sigma_t)
            s[j,i,2]=rnorm(1,mu[j,2],sigma_t)
            if(s[j,i,1] < xlim[2]+buff & s[j,i,1] > xlim[1]-buff & s[j,i,2] < ylim[2]+buff & s[j,i,2] > ylim[1]-buff){
              out=0
            }
            countout=countout+1
            if(countout==5000){#So you don't end up in an infinite loop if you set sigma_t too large relative to state space size
              stop("5000 ACs proposed outside of state space. Should probably reconsider settings.")
            }
          }
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }
    }
    psi=N[1]/M
    #Calculate (per capita) gamma or N
    if(length(phi)==1){
      phi=rep(phi,t-1)
    }

    if(is.null(gamma)){
      gamma=rep(NA,t-1)
      for(l in 2:t){
        gamma[l-1]=(N[l]-phi[t-1]*N[l-1])/(N[l-1])
      }
      storeparms$gamma=gamma
    }else{
      if(length(gamma)==1){
        gamma=rep(gamma,t-1)
      }
      N=c(N,rep(NA,t-1))
      for(l in 2:t){
        N[l]=round(N[l-1]*phi[t-1]+N[l-1]*gamma[l-1])
      }
    }
    #####Population Dynamics############# Stolen from Richard Chandler
    z=a=matrix(NA,M,t)
    # z[,1] <- rbinom(M, 1, EN1/M)
    z[1:N[1],1]=1
    z[(N[1]+1):M]=0
    a[,1]= 1-z[,1] # Available to be recruited?
    gamma.prime=rep(NA,4)
    for(i in 2:t) {
      ER <- sum(z[,i-1])*gamma[i-1] # Expected number of recruits
      A <- sum(a[,i-1]) # nAvailable to be recruited
      if(ER>A){
        stop("M is too low. There aren't any individuals left to be recruited")
      }
      gamma.prime[i-1] <- ER/A # individual-level recruitment *probability*
      Ez <- z[,i-1]*phi[i-1] + a[,i-1]*gamma.prime[i-1]
      z[,i] <- rbinom(M, 1, Ez)
      a[,i] <- apply(z[,1:i]<1, 1, all)
    }

    #######Capture process######################
    # Simulate encounter history
    idx=0
    y <-array(0,dim=c(M,maxJ,t))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            y[i,j,l]=rbinom(1,K[l],pd[i,j,l]*z[i,l])
          }
        }
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            y[i,j,l]=rpois(1,K[l]*lamd[i,j,l]*z[i,l])
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    caps=apply(y,1,sum)
    idx=order(caps,decreasing=TRUE)
    y=y[idx,,]
    s=s[idx,,]
    mu=mu[idx,]
    keep=which(rowSums(y)>0)
    y=y[keep,,]
    s=s[keep,,]
    mu=mu[keep,]
    n=sum(caps>0)
    caps2d=apply(y,c(1,3),sum)
    n2d=colSums(caps2d>0)
    if(length(storeparms$gamma)==1){
      gamma=gamma[1]
    }else{
      gamma=gamma
    }
    if(is.null(sigma_t)){
      out<-list(y=y,s=s,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype)
    }else{
      out<-list(y=y,mu=mu,s=s,X=X,K=K,n=n,n2d=n2d,buff=buff,J=J,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype)
    }
    return(out)
  }
