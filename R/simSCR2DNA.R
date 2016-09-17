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


#' Simulate data from a SCR study with DNA collected from scat and hair snares (update help later)
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

simSCR2DNA <-
  function(N=120,lam01=0.2,lam01b=0,lam02=0.2,sigma=0.50,K1=10,K2=10,X1=X1,X2=X2,buff=2){
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(c(X1[,1],X2[,1]))-buff,max(c(X1[,1],X2[,1]))+buff), runif(N,min(c(X1[,2]),X2[,2])-buff,max(c(X1[,2]),X2[,2])+buff))
    D1<- e2dist(s,X1)
    D2<- e2dist(s,X2)
    lamd1<- lam01*exp(-D1*D1/(2*sigma*sigma))
    lamd2<- lam02*exp(-D2*D2/(2*sigma*sigma))
    J1<- nrow(X1)
    J2<- nrow(X2)
    # Simulate hair encounter history
    y1 <-array(0,dim=c(N,K1,J1))
    pd1=cellprobsSCR(lamd1)
    if(lam01b==0){ #if no behavioral response
      for(i in 1:N){
        for(j in 1:J1){
          for(k in 1:K1){
            y1[i,k,j]=rbinom(1,1,pd1[i,j])
          }
        }
      }
    }else{
      lamd1b=lam01b*exp(-D1*D1/(2*sigma*sigma))
      pd1b=cellprobsSCR(lamd1b)
      state=matrix(0,nrow=N,ncol=J1) #Matrix of indices 1 indicating previously captured at trap 0 o.w.
      for(i in 1:N){
        for(j in 1:J1){
          for(k in 1:K1){
            y1[i,k,j]=rbinom(1,1,pd1[i,j]*(1-state[i,j])+pd1b[i,j]*state[i,j])
            state[i,j]=max(state[i,j],y1[i,k,j])  #update state
          }
        }
      }
    }
    # Simulate scat encounter history
    y2 <-array(0,dim=c(N,K2,J2))
    pd2=cellprobsSCR(lamd2)
    for(i in 1:N){
      for(j in 1:J2){
        for(k in 1:K2){
          y2[i,k,j]=rbinom(1,1,pd2[i,j])
        }
      }
    }
    #Remove guys not captured.  or not for now
    rem=which(rowSums(y1)==0&rowSums(y2)==0)
    # y1=y1[-rem,,]
    # y2=y2[-rem,,]
    n=c(nrow(y1)-length(rem),length(which(rowSums(y1)>0)),length(which(rowSums(y2)>0)))
    out<-list(y1=y1,y2=y2,s=s,X1=X1,X2=X2,K1=K1,K2=K2,n=n,buff=buff)
    return(out)
  }
