#' Calculate goodness of fit statistics for a model fit by mcmc.genCatSMR
#' @param data a list produced by simSCR2DNA or in the same format
#' @param posterior a posterior produced by mcmc.SCR
#' @param  use a vector if iteration values at which the statistics should be calculated
#' @author Ben Augustine
#' @description Calculate goodness of fit statistics for a model fit by mcmc.SCR
#' @export

# posterior=out
# use=seq(1,5333,10)
# obstype=c("bernoulli","poisson")
genCatSMR.PPcheck=function(data,posterior,use,obstype,tf1=NA,tf2=NA){
  X1=data$X1
  X2=data$X2
  J1=dim(data$y.mark)[2]
  J2=dim(data$y.sight.marked)[2]
  K1=dim(data$y.mark)[3]
  K2=dim(data$y.sight.marked)[3]
  if(any(is.na(tf1))){
    K2D1=matrix(rep(tf1,M),nrow=M,ncol=J1,byrow=TRUE)
  }else{
    K2D1=tf1
  }
  if(any(is.na(tf2))){
    K2D2=matrix(rep(tf1,M),nrow=M,ncol=J2,byrow=TRUE)
  }else{
    K2D2=tf2
  }
  
  M=dim(posterior$zout)[2]
  n.marked=dim(data$y.mark)[1]
  n.samples=nrow(data$y.sight.unmarked)
  y.mark.ijobs=apply(data$y.mark,c(1,2),sum)
  y.mark.ijobs=rbind(y.mark.ijobs,matrix(0,nrow=M-n.marked,ncol=J1))
  y.mark.iobs=rowSums(y.mark.ijobs)
  y.mark.jobs=colSums(y.mark.ijobs)
  #statistic structures
  niter=length(use)
  n.samples.prime=n.um.prime=rep(NA,niter)
  T1.mark.obs=T1.mark=rep(NA,niter)
  T2.mark.obs=T2.mark=matrix(NA,nrow=niter,ncol=M)
  T3.mark.obs=T3.mark=matrix(NA,nrow=niter,ncol=J1)
  
  # y.sight.ijobs=apply(data$y,c(1,2),sum)
  # y.sight.jobs=rbind(y.sight.ijobs,matrix(0,nrow=M-n.marked,ncol=J2))
  # y.sight.ijobs=rbind(y.sight.ijobs,matrix(0,nrow=M-n.marked,ncol=J2))
  # y.sight.iobs=rowSums(y.sight.ijobs)
  # y.sight.jobs=colSums(y.sight.ijobs)
  T1.sight.obs=T1.sight=rep(NA,niter)
  T2.sight.obs=T2.sight=matrix(NA,nrow=niter,ncol=M)
  T3.sight.obs=T3.sight=matrix(NA,nrow=niter,ncol=J2)
  # T2.sight.UM=T2.sight.obs.UM=rep(NA,niter)

  idx=1
  for(iter in use){
    #extract parameter values
    lam0.mark=posterior$out[iter,1]
    lam0.sight=posterior$out[iter,2]
    sigma=posterior$out[iter,3]
    if(dim(posterior$out)[2]==6){
      s1=s2=cbind(posterior$sxout[iter,],posterior$syout[iter,])
    }else{
      s1=cbind(posterior$s1xout[iter,],posterior$s1yout[iter,])
      s2=cbind(posterior$s2xout[iter,],posterior$s2yout[iter,])
    }
    z=posterior$zout[iter,]
    ID=posterior$IDout[iter,]
    
    #reconstruct y.sight
    y.sight=array(0,dim=c(M,J2,K2))
    y.sight[1:n.marked,,]=data$y.sight.marked
    for(i in 1:n.samples){
      y.sight[ID[i],,]=y.sight[ID[i],,]+data$y.sight.unmarked[i,,]
    }
    y.sight.ijobs=apply(y.sight,c(1,2),sum)
    y.sight.iobs=rowSums(y.sight.ijobs)
    y.sight.jobs=colSums(y.sight.ijobs)
    
    #Marking process
    D1<- e2dist(s1,X1)
    lamd.mark<- lam0.mark*exp(-D1^2/(2*sigma*sigma))
    y.mark <-array(0,dim=c(M,J1))
    if(obstype[1]=="bernoulli"){
      pd.mark=1-exp(-lamd.mark)
      for(i in 1:M){
        for(j in 1:J1){
            y.mark[i,j]=rbinom(1,K2D1[i,j],z[i]*pd.mark[i,j]) 
          }
      }
    }else if(obstype[1]=="poisson"){
      for(i in 1:M){
        for(j in 1:J1){
            y.mark[i,j]=rpois(1,K2D1[i,j]*z[i]*lamd.mark[i,j])
          }
      }
    }else{
      stop("obstype 1 not recognized")
    }
    #Sighting process
    D2<- e2dist(s2,X2)
    lamd.sight<- lam0.sight*exp(-D2^2/(2*sigma*sigma))
    y.sight <-array(0,dim=c(M,J2))
    if(obstype[2]=="bernoulli"){
      pd.sight=1-exp(-lamd.sight)
      for(i in 1:M){
        for(j in 1:J2){
            y.sight[i,j]=rbinom(1,K2D2[i,j],z[i]*pd.sight[i,j]) 
        }
      }
    }else if(obstype[2]=="poisson"){
      for(i in 1:M){
        for(j in 1:J2){
          y.sight[i,j]=rpois(1,z[i]*K2D2[i,j]*lamd.sight[i,j])
        }
      } 
    }else{
      stop("obstype 2 not recognized")
    }
    ##marking process stats
    #T1 - ind by trap probs
    if(obstype[1]=="bernoulli"){
      Eyij=pd.mark*K2D1
    }else{
      Eyij=lamd.mark*K2D1
    }
    T1.mark[idx]=sum((sqrt(y.mark)-sqrt(Eyij))^2)
    T1.mark.obs[idx]=sum((sqrt(y.mark.ijobs)-sqrt(Eyij))^2)
    #T2 - individual probs
    Eyi=rowSums(Eyij)
    yi=rowSums(y.mark)
    # T2.mark[idx]=sum((sqrt(yi)-sqrt(Eyi))^2)
    # T2.mark.obs[idx]=sum((sqrt(y.mark.iobs)-sqrt(Eyi))^2)
    T2.mark[idx,]=((sqrt(yi)-sqrt(Eyi))^2)
    T2.mark.obs[idx,]=((sqrt(y.mark.iobs)-sqrt(Eyi))^2)
    #T3 - trap probs
    Eyj=colSums(Eyij)
    yj=colSums(y.mark)
    # T3.mark[idx]=sum((sqrt(yj)-sqrt(Eyj))^2)
    # T3.mark.obs[idx]=sum((sqrt(y.mark.jobs)-sqrt(Eyj))^2)
    T3.mark[idx,]=((sqrt(yj)-sqrt(Eyj))^2)
    T3.mark.obs[idx,]=((sqrt(y.mark.jobs)-sqrt(Eyj))^2)
    
    ##sighting process stats
    #T1 - ind by trap probs
    if(obstype[2]=="bernoulli"){
      Eyij=pd.sight*K2D2
    }else{
      Eyij=lamd.sight*K2D2
    }
    T1.sight[idx]=sum((sqrt(y.sight)-sqrt(Eyij))^2)
    T1.sight.obs[idx]=sum((sqrt(y.sight.ijobs)-sqrt(Eyij))^2)
    #T2 - individual probs
    Eyi=rowSums(Eyij)
    yi=rowSums(y.sight)
    # T2.sight[idx]=sum((sqrt(yi)-sqrt(Eyi))^2)
    # T2.sight.obs[idx]=sum((sqrt(y.sight.iobs)-sqrt(Eyi))^2)
    T2.sight[idx,]=((sqrt(yi)-sqrt(Eyi))^2)
    T2.sight.obs[idx,]=((sqrt(y.sight.iobs)-sqrt(Eyi))^2)
    # T2.sight.UM[idx]=sum(T2.sight[idx,unique(ID)])
    # T2.sight.obs.UM[idx]=sum(T2.sight.obs[idx,unique(ID)])
    #T3 - trap probs
    Eyj=colSums(Eyij)
    yj=colSums(y.sight)
    # T3.sight[idx]=sum((sqrt(yj)-sqrt(Eyj))^2)
    # T3.sight.obs[idx]=sum((sqrt(y.sight.jobs)-sqrt(Eyj))^2)
    T3.sight[idx,]=((sqrt(yj)-sqrt(Eyj))^2)
    T3.sight.obs[idx,]=((sqrt(y.sight.jobs)-sqrt(Eyj))^2)
    
    n.samples.prime[idx]=sum(y.sight[n.marked:nrow(y.sight),])
    n.um.prime[idx]=sum(rowSums(y.sight[n.marked:nrow(y.sight),])>0)
    idx=idx+1
  }
  n.samples.per.um.prime=n.samples.prime/n.um.prime
  n.samples.obs=nrow(data$y.sight.unmarked)
  if(dim(posterior$out)[2]==6){
    n.um.obs=out$out[use,5]
  }else{
    n.um.obs=out$out[use,6]
  }
  n.samples.per.um.obs=n.samples.obs/n.um.obs
  allstats=list(n.samples.prime,n.um.prime,n.samples.per.um.prime,
                      n.samples.obs,n.um.obs,n.samples.per.um.obs,T1.mark,T2.mark,T3.mark,
                      T1.sight,T2.sight,T3.sight,T1.mark.obs,T2.mark.obs,T3.mark.obs,
                      T1.sight.obs,T2.sight.obs,T3.sight.obs)
  names(allstats)=c("n.samples.prime","n.um.prime","n.samples.per.um.prime",
                    "n.samples.obs","n.um.obs","n.samples.per.um.obs","T1.mark","T2.mark","T3.mark",
                    "T1.sight","T2.sight","T3.sight","T1.mark.obs","T2.mark.obs","T3.mark.obs",
                    "T1.sight.obs","T2.sight.obs","T3.sight.obs")
  
  return(allstats)
}
