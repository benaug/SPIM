#' Simulate data from a SCR study and simulate a random trap failure process
#' @param N a vector indicating the number of individuals to simulate
#' @param lam0 the detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the K x 2 matrix of trap locations
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param failprob the probability a trap will fail on each occasion
#' @param faildur the number of occasions a trap remains inoperable
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description See the simSCR help file for a description of the capture process.  The trap failure process operates
#' like this:  traps fail with probabilty failprob and remain disabled for faildur occasions. This code does not keep up with
#' multiple cameras at a site like the 2side simulation code.  A trap is either on or off.
#' @author Ben Augustine
#' @export

simSCRtf <-
  function(N=120,lam0=0.2,lam0b=0,sigma=0.70,K=10,X=X,buff=3,failrate=0.05,faildur=2,obstype="bernoulli"){
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    #Simulate trap failure/correction process
    onoff=matrix(NA,nrow=J,ncol=K)
    onoff[,1]=1
    offoccs=rep(0,J)
    on=rep(1,J)
    for(k in 2:K){
      offoccs[offoccs>0]=offoccs[offoccs>0]+1 #if off, add another occasion
      onoff[,k]=rbinom(J,on,1-failrate)
      offoccs[(onoff[,k]==0)&(on==1)]=1 #Turn off stations that just failed
      on[(onoff[,k]==0)&(on==1)]=0 #Turn off stations that just failed
      backon=which(offoccs==faildur)
      if(length(backon)>0){
        on[backon]=1
        offoccs[backon]=0
      }
    }
    # Simulate encounter history
    y <-array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=cellprobsSCR(lamd)
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,j,k]=rbinom(1,1,pd[i,j]*onoff[j,k])
            }
          }
        }
    }else if(obstype=="poisson"){
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,j,k]=rbinom(1,1,lamd[i,j]*onoff[j,k])
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
    tf=rowSums(onoff)
    out<-list(y=y,s=s,X=X, K=K,n=n,buff=buff,tf=tf,obstype=obstype)
    return(out)
  }
