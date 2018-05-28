#' Calculate goodness of fit statistics for a model fit by SCR2DNAmcmc
#' @param data a list produced by simSCR2DNA or in the same format
#' @param posterior a posterior produced by SCR2DNAmcmc
#' @param  use a vector if iteration values at which the statistics should be calculated
#' @param sharesig a logical indicating whether or not sigma was shared
#' @author Ben Augustine
#' @description Calculate goodness of fit statistics for a model fit by SCR2DNAmcmc
#' @export
SCR2DNAPPcheck=function(data,posterior,use,sharesig){
  X1=data$X1
  X2=data$X2
  K1=dim(data$y1)[3]
  K2=dim(data$y2)[3]
  J1=dim(data$y1)[2]
  J2=dim(data$y2)[2]
  M=dim(posterior$zout)[2]
  n=data$n
  
  y1ijobs=apply(data$y1,c(1,2),sum)
  y2ijobs=apply(data$y2,c(1,2),sum)
  y1ijobs=rbind(y1ijobs,matrix(0,nrow=M-n,ncol=J1))
  y2ijobs=rbind(y2ijobs,matrix(0,nrow=M-n,ncol=J2))
  y1iobs=rowSums(y1ijobs)
  y2iobs=rowSums(y2ijobs)
  y1jobs=colSums(y1ijobs)
  y2jobs=colSums(y2ijobs)
  
  
  niter=length(use)
  T11obs=T12obs=T11=T12=rep(NA,niter)
  T21obs=T22obs=T21=T22=rep(NA,niter)
  T31obs=T32obs=T31=T32=rep(NA,niter)
  ncap1=ncap2=ncap3=bothcap=bothcap2=rep(NA,niter)
  T21full=T21fullobs=matrix(NA,nrow=M,ncol=niter)
  T22full=T22fullobs=matrix(NA,nrow=M,ncol=niter)
  T31full=T31fullobs=matrix(NA,nrow=127,ncol=niter)
  T32full=T32fullobs=matrix(NA,nrow=98,ncol=niter)

  
  idx=1
  for(iter in use){
    lam01=posterior$out[iter,1]
    lam02=posterior$out[iter,2]
    s=cbind(posterior$sxout[iter,],posterior$syout[iter,])
    z=posterior$zout[iter,]
    if(sharesig){
      sigma=rep(posterior$out[iter,3],2)
    }else{
      sigma=c(posterior$out[iter,3],posterior$out[iter,4])
    }
    
    #######Capture process######################
    D1<- e2dist(s,X1)
    D2<- e2dist(s,X2)
    lamd1<- lam01*exp(-D1*D1/(2*sigma[1]*sigma[1]))
    lamd2<- lam02*exp(-D2*D2/(2*sigma[2]*sigma[2]))
    pd1=1-exp(-lamd1)
    pd2=1-exp(-lamd2)
    J1<- nrow(X1)
    J2<- nrow(X2)
    # Simulate hair encounter history
    y1 <-array(0,dim=c(M,K1,J1))
    for(i in 1:M){
      for(j in 1:J1){
        for(k in 1:K1){
          y1[i,k,j]=rbinom(1,1,z[i]*pd1[i,j])
        }
      }
    }
    # Simulate scat encounter history
    y2 <-array(0,dim=c(M,K2,J2))
    for(i in 1:M){
      for(j in 1:J2){
        for(k in 1:K2){
          y2[i,k,j]=rbinom(1,1,z[i]*pd2[i,j])
        }
      }
    }
    #T1 - ind by trap probs
    Eyij1=pd1*K1
    Eyij2=pd2*K2
    y1ij=apply(y1,c(1,3),sum)
    y2ij=apply(y2,c(1,3),sum)
    T11[idx]=sum((sqrt(y1ij)-sqrt(Eyij1))^2)
    T12[idx]=sum((sqrt(y2ij)-sqrt(Eyij2))^2)
    T11obs[idx]=sum((sqrt(y1ijobs)-sqrt(Eyij1))^2)
    T12obs[idx]=sum((sqrt(y2ijobs)-sqrt(Eyij2))^2)
    #T2 - individual probs
    Eyi1=rowSums(Eyij1)
    Eyi2=rowSums(Eyij2)
    y1i=rowSums(y1ij)
    y2i=rowSums(y2ij)
    T21[idx]=sum((sqrt(y1i)-sqrt(Eyi1))^2)
    T22[idx]=sum((sqrt(y2i)-sqrt(Eyi2))^2)
    T21obs[idx]=sum((sqrt(y1iobs)-sqrt(Eyi1))^2)
    T22obs[idx]=sum((sqrt(y2iobs)-sqrt(Eyi2))^2)
    #T2 full
    T21full[,idx]=((sqrt(y1i)-sqrt(Eyi1))^2)
    T22full[,idx]=((sqrt(y2i)-sqrt(Eyi2))^2)
    T21fullobs[,idx]=((sqrt(y1iobs)-sqrt(Eyi1))^2)
    T22fullobs[,idx]=((sqrt(y2iobs)-sqrt(Eyi2))^2)
    #T3 - trap probs
    Eyj1=colSums(Eyij1)
    Eyj2=colSums(Eyij2)
    y1j=colSums(y1ij)
    y2j=colSums(y2ij)
    T31[idx]=sum((sqrt(y1j)-sqrt(Eyj1))^2)
    T32[idx]=sum((sqrt(y2j)-sqrt(Eyj2))^2)
    T31obs[idx]=sum((sqrt(y1jobs)-sqrt(Eyj1))^2)
    T32obs[idx]=sum((sqrt(y2jobs)-sqrt(Eyj2))^2)
    #T3 full
    T31full[,idx]=((sqrt(y1j)-sqrt(Eyj1))^2)
    T32full[,idx]=((sqrt(y2j)-sqrt(Eyj2))^2)
    T31fullobs[,idx]=((sqrt(y1jobs)-sqrt(Eyj1))^2)
    T32fullobs[,idx]=((sqrt(y2jobs)-sqrt(Eyj2))^2)
    #ncaps
    ncap1[idx]=sum(rowSums(y1)>0)
    ncap2[idx]=sum(rowSums(y2)>0)
    ncap3[idx]=sum(rowSums(y1)>0&rowSums(y2)>0)
    idx=idx+1
  }

  allstats=data.frame(T11=T11,T12=T12,T21=T21,T22=T22,T31=T31,T32=T32,
                      T11obs=T11obs,T12obs=T12obs,T21obs=T21obs,T22obs=T22obs,T31obs=T31obs,T32obs=T32obs,
                      ncap1=ncap1,ncap2=ncap2,ncap3=ncap3)
  return(allstats)
  
}