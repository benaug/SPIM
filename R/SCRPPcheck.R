#' Calculate goodness of fit statistics for a model fit by mcmc.SCR
#' @param data a list produced by simSCR2DNA or in the same format
#' @param posterior a posterior produced by mcmc.SCR
#' @param  use a vector if iteration values at which the statistics should be calculated
#' @author Ben Augustine
#' @description Calculate goodness of fit statistics for a model fit by mcmc.SCR
#' @export
SCRPPcheck=function(data,posterior,use){
  X=data$X
  K=dim(data$y)[2]
  J=dim(data$y)[3]
  M=dim(posterior$zout)[2]
  n=data$n
  rem=which(rowSums(data$y)==0)
  if(length(rem)>0){
    data$y=data$y[-rem,,]
  }
  yijobs=apply(data$y,c(1,3),sum)
  yjobs=rbind(yijobs,matrix(0,nrow=M-n,ncol=J))
  yijobs=rbind(yijobs,matrix(0,nrow=M-n,ncol=J))
  yiobs=rowSums(yijobs)
  yjobs=colSums(yijobs)
  niter=length(use)
  T1obs=T1=rep(NA,niter)
  T2obs=T2=rep(NA,niter)
  T3obs=T3=rep(NA,niter)
  idx=1
  for(iter in use){
    #extract parameter values
    lam0=posterior$out[iter,1]
    sigma=rep(posterior$out[iter,2],2)
    s=cbind(posterior$sxout[iter,],posterior$syout[iter,])
    z=posterior$zout[iter,]
    
    # Simulate encounter history
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    pd=1-exp(-lamd)
    y <-array(0,dim=c(M,K,J))
    for(i in 1:M){
      for(j in 1:J){
        for(k in 1:K){
          y[i,k,j]=rbinom(1,1,z[i]*pd[i,j])
        }
      }
    }
    #T1 - ind by trap probs
    Eyij=pd*K
    yij=apply(y,c(1,3),sum)
    T1[idx]=sum((sqrt(yij)-sqrt(Eyij))^2)
    T1obs[idx]=sum((sqrt(yijobs)-sqrt(Eyij))^2)
    #T2 - individual probs
    Eyi=rowSums(Eyij)
    yi=rowSums(yij)
    T2[idx]=sum((sqrt(yi)-sqrt(Eyi))^2)
    T2obs[idx]=sum((sqrt(yiobs)-sqrt(Eyi))^2)
    #T3 - trap probs
    Eyj=colSums(Eyij)
    yj=colSums(yij)
    T3[idx]=sum((sqrt(yj)-sqrt(Eyj))^2)
    T3obs[idx]=sum((sqrt(yjobs)-sqrt(Eyj))^2)
    idx=idx+1
  }
  allstats=data.frame(T1=T1,T2=T2,T3=T3,T1obs=T1obs,T2obs=T2obs,T3obs=T3obs)
  return(allstats)
}
