#' Simulate data from a SCR study with DNA collected from scat and hair snares
#' @param N a vector indicating the number of individuals to simulate
#' @param lam01 the hair snare detection function hazard rate
#' @param lam0b the hair snare detection function hazard rate after capture at a trap
#' @param lam02 the scat detection function hazard rate
#' @param sigma a vector of spatial scale parameters for the hair snare and scat detection functions (in that order)
#' @param K1 the number of hair snare capture occasions
#' @param K2 the number of scat capture occasions
#' @param X1 the K1 x 2 matrix of hair snare locations
#' @param X2 the K2 x 2 matrix of scat collection locations (descretized transects)
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @return a list containing the capture histories (y1 for hair snares, y2 for scat collection), activity centers, trap objects, and several other data objects and summaries.
#' @description This function simulates data from a joint hair snare and scat collection SCR study. The extent of the state space is controlled by "buff", which buffers the
#' minimum and maximum X and Y extents.  Therefore, it probably only makes sense for square or rectangular grids.  Functionality
#' for user-provided polygonal state spaces will be added in the future.
#' @author Ben Augustine
#' @export

simSCR2DNA <-
  function(N=120,lam01=0.2,lam01b=0,lam02=0.2,sigma=c(0.5,0.5),K1=10,K2=10,X1=X1,X2=X2,buff=2){
    if(length(sigma)!=2){
      stop("sigma vector must be of length 2")
    }
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(c(X1[,1],X2[,1]))-buff,max(c(X1[,1],X2[,1]))+buff), runif(N,min(c(X1[,2]),X2[,2])-buff,max(c(X1[,2]),X2[,2])+buff))
    D1<- e2dist(s,X1)
    D2<- e2dist(s,X2)
    lamd1<- lam01*exp(-D1*D1/(2*sigma[1]*sigma[1]))
    lamd2<- lam02*exp(-D2*D2/(2*sigma[2]*sigma[2]))
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
      lamd1b=lam01b*exp(-D1*D1/(2*sigma[1]*sigma[1]))
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
