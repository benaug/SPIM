SCRmcmcRcpp <-
  function(data,niter=1000,nburn=0, nthin=1, M =NA,K=NA, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           storeLatent=TRUE,obstype="bernoulli"){
    library(abind)
    y<-data$y
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- data$K
    n=dim(y)[1]
    #Reduce to 2D data
    if(length(dim(y))==3){
      y=apply(y,c(1,2),sum)
    }
    
    #If using polygon state space
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
      xlim=c(min(vertices[,1]),max(vertices[,1]))
      ylim=c(min(vertices[,2]),max(vertices[,2]))
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
      ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
      vertices=cbind(xlim,ylim)
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    #Trap operation
    if("tf"%in%names(data)){
      tf=data$tf
      if(length(tf)!=nrow(X)){
        stop("If using a trap operation file, must input vector of length nrow(X).")
      }
      
    }else{
      tf=rep(K,J)
    }
    #Don't do this for Rcpp
    # tf=matrix(rep(tf,M),ncol=J,nrow=M,byrow=TRUE)
    
    ##pull out initial values
    psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    known.vector=c(rep(1,n),rep(0,M-n))
    y=abind(y,matrix(0,nrow=M-n,ncol=J),along=1)
    #Initialize z
    z=1*(apply(y,1,sum)>0)
    add=M*(0.5-sum(z)/M)
    if(add>0){
      z[sample(which(z==0),add)]=1 #switch some uncaptured z's to 1.
    }
    
    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y[i,]>0,1:2],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s[i,]<- trps
      }
    }
    if(useverts==TRUE){
      inside=rep(NA,nrow(s))
      for(i in 1:nrow(s)){
        inside[i]=inout(s[i,],vertices)
      }
      idx=which(inside==FALSE)
      if(length(idx)>0){
        for(i in 1:length(idx)){
          while(inside[idx[i]]==FALSE){
            s[idx[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside[idx[i]]=inout(s[idx[i],],vertices)
          }
        }
      }
    }
    D=e2dist(s, X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y=array(0,dim=c(M,J))
    #check if ll finite
    tf2D=matrix(rep(tf,M),ncol=J,nrow=M,byrow=TRUE)
    
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      ll.y=dbinom(y,tf2D,pd*z,log=TRUE)
    }else if(obstype=="poisson"){
      ll.y=dpois(y,tf2D*lamd*z,log=TRUE)
    }
    ll.y.cand=ll.y
    
    # some objects to hold the MCMC output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","psi"))
    if(storeLatent){
      sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
    }
    if(obstype=="bernoulli"){
      obstypein=1
    }else{
      obstypein=2
    }
    if(!is.finite(sum(ll.y)))stop("Observation model likelihood is not finite. Try raising starting values for lam0 and/or sigma.")
    
  
#Run MCMC
store=MCMC1(lam0,sigma,y,lamd,z,X,K,D,known.vector,s,
                  psi,xlim,ylim,useverts,vertices,
                  proppars$lam0,proppars$sigma,proppars$sx,
                      proppars$sy,niter,nburn,nthin,obstypein,
                  tf,storeLatent)

if(storeLatent){
  list(out=store[[1]], sxout=store[[2]], syout=store[[3]], zout=store[[4]],data=data)
}else{
  list(out=store[[1]],data=data)
}
}
