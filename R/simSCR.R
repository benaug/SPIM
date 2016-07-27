e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

cellprobsSCR<- function(lamd){
  # For gaussian hazard model convert lamda(s,x) to p(s,x)
  N<- dim(lamd)[1]
  J<- dim(lamd)[2] # traps
  pmat<- matrix(NA,nrow=N,ncol=J)
  for(j in 1:J){
    pmat[,j]<- 1-exp(-lamd[,j])
  }
  pmat
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

simSCR <-
  function(N=120,lam0=0.2,lam0b=0,sigma=0.50,K=10,X=X,buff=3,obstype="bernoulli"){
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    # Simulate encounter history
    y <-array(0,dim=c(N,K,J))
    if(obstype=="bernoulli"){
      pd=cellprobsSCR(lamd)
      if(lam0b==0){ #if no behavioral response
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,pd[i,j])
            }
          }
        }
      }else{
        if((-lam0b)>lam0){stop("-b must be > lam0")}
        lamdb=(lam0+lam0b)*exp(-D*D/(2*sigma*sigma))
        pdb=cellprobsSCR(lamdb)
        state=matrix(0,nrow=N,ncol=J) #Matrix of indices 1 indicating previously captured at trap 0 o.w.
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,pd[i,j]*(1-state[i,j])+pdb[i,j]*state[i,j])
              state[i,j]=max(state[i,j],y[i,k,j])  #update state
            }
          }
        }
      }
    }else if(obstype=="poisson"){
      if(lam0b==0){ #if no behavioral response
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,lamd[i,j])
            }
          }
        }
      }else{
        if((-lam0b)>lam0){stop("-b must be > lam0")}
        lamdb=(lam0+lam0b)*exp(-D*D/(2*sigma*sigma))
        state=matrix(0,nrow=N,ncol=J) #Matrix of indices 1 indicating previously captured at trap 0 o.w.
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rpois(1,lamd[i,j]*(1-state[i,j])+lamdb[i,j]*state[i,j])
              state[i,j]=max(state[i,j],y[i,k,j])  #update cap
            }
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    caps=apply(y,1,sum)
    idx=order(caps,decreasing=TRUE)
    y=y[idx,,]
    s=s[idx,]
    n=sum(caps>0)
    y=y[rowSums(y)>0,,]
    #Count spatial recaps
    y2D=apply(y,c(1,3),sum)
    scaps=rowSums(1*(y2D>0))
    scaps[scaps>0]=scaps[scaps>0]-1
    nscap=sum(scaps>0)
    sumscap=sum(scaps)
    out<-list(y=y,s=s,X=X, K=K,n=n,nscap=nscap,sumscap=sumscap,buff=buff,obstype=obstype)
    return(out)
  }
