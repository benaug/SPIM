# Run MCMC algorithm for SCR model with inhomogenous density and 1 categorical density covariate. Discrete state space. In Rcpp.
# @param data a list produced by simSCRIPP or in the same format
# @param niter number of MCMC iterations to run
# @param  nburn number of MCMC iterations to discard as burn in
# @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
# @param M The size of the augmented superpopulation
# @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam01=lam01,lam02=lam02,sigma=sigma)
# @param proppars a list of tuning parameters for the proposal distributions
# @return  a list with the posteriors for the SCR detection function and density covariate parameters (out), s, z
# @description code not yet functional
# @author Ben Augustine
# @export

SCRIPPmcmcRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,s=1,beta0=0.05,beta1=0.05),cellArea,keepACs=FALSE){
    ###
    library(abind)
    y<-data$y
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(y)[2]
    n<- data$n
    grid=as.matrix(data$grid)
    npix=nrow(grid)
    ##pull out initial values
    beta0<- inits$beta0
    beta1<- inits$beta1
    lam0<- inits$lam0
    sigma<- inits$sigma

    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],K, J)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))

    #Make initial complete data set
    y2D=apply(y,c(1,3),sum)

    #Optimize starting locations given where they are trapped.
    cell=sample(1:npix,M)#assign random locations
    s<- cbind(cell,grid[cell,1:2])
    idx=which(rowSums(y)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[y2D[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      loc=c(mean(trps[,1]),mean(trps[,2]))
      d=sqrt((loc[1]-grid[,1])^2+(loc[2]-grid[,2])^2)
      cell=which(d==min(d))
      if(length(cell)>1){
        cell=cell[1]
      }
      s[i,]=c(cell,grid[cell,1:2])
    }

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=6)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","EN","beta0","beta1"))
    sxout<- syout<- scellout<-zout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    D<- e2dist(s[,2:3], X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))

    mu=exp(beta0+beta1*1*(grid[,3]==2))*cellArea
    EN=sum(mu)
    psi <- EN/M
    if (psi > 1){
      stop("Bad initial values for beta0 or beta1. Or M is too low")
    }
    z=1*(apply(y2D,1,sum)>0)
    N=sum(z)
    N2=round(psi*EN)

    if(N<N2){
      z[sample(which(z==0),N2-N)]=1 #switch some uncaptured z's to 1.
      N=sum(z)
    }

    #Run MCMC
    store=MCMC1b(lam0,sigma,beta0,beta1,y2D,z,X,K,D,N,known.vector,s[,2:3],s[,1],psi,grid,cellArea,EN,proppars$lam0,proppars$sigma,
                 proppars$beta0,proppars$beta1,proppars$sx,proppars$sy,niter,nburn,nthin)


    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]],zout=store[[4]])
    }else{
      list(out=store[[1]])
    }
  }
