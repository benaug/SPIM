SCRmcmcOpenRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    library(abind)
    t=dim(data$y)[3]
    y<-data$y
    X<-data$X
    J<-data$J
    maxJ=max(J)
    K<-data$K
    n=dim(data$y)[1]
    n2d<-colSums(apply(data$y,c(3),rowSums)>0)
    metamu="sigma_t"%in%names(inits)
    ####Error checks
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    #If using polygon state space
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
      xlim=c(min(vertices[,1]),max(vertices[,1]))
      ylim=c(min(vertices[,2]),max(vertices[,2]))
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(unlist(lapply(X,function(x){min(x[,1])}))),min(unlist(lapply(X,function(x){max(x[,1])}))))+c(-buff, buff)
      ylim<- c(min(unlist(lapply(X,function(x){min(x[,2])}))),min(unlist(lapply(X,function(x){max(x[,2])}))))+c(-buff, buff)
      vertices=rbind(c(xlim[1],ylim[1]),c(xlim[1],ylim[2]),c(xlim[2],ylim[2]),c(xlim[2],ylim[1]))
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    ##pull out initial values
    # psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    sigma_t=inits$sigma_t
    gamma=inits$gamma
    phi=inits$phi
    psi=inits$psi
    if(!length(lam0)%in%c(1,t)){
      stop("Input either 1 or t initial values for lam0")
    }
    if(!length(sigma)%in%c(1,t)){
      stop("Input either 1 or t initial values for sigma")
    }
    if(!length(gamma)%in%c(1,t-1)){
      stop("Input either 1 or t-1 initial values for gamma (psi and fixed gamma)")
    }
    if(!length(phi)%in%c(1,t-1)){
      stop("Input either 1 or t-1 initial values for phi")
    }
    #Check proppars
    if(length(proppars$propz)!=(t-1)){
      stop("must supply t-1 proppars for propz")
    }
    if(length(lam0)!=length(proppars$lam0)){
      stop("Must supply a tuning parameter for each lam0")
    }
    if(length(sigma)!=length(proppars$sigma)){
      stop("Must supply a tuning parameter for each sigma")
    }
    if(length(gamma)!=length(proppars$gamma)){
      stop("Must supply a tuning parameter for each gamma")
    }
    if(length(phi)!=length(proppars$phi)){
      stop("Must supply a tuning parameter for each phi")
    }
    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],maxJ, t)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))

    #Initialize z, r, and a consistent with y
    known.matrix=1*(apply(y,c(1,3),sum)>0)#make z consistent with y
    if(t>2){
      for(l in 2:(t-1)){#Turn on zeros with 1's on either side
        known.matrix[known.matrix[,l]==0&known.matrix[,l-1]==1&rowSums(matrix(known.matrix[,(l+1):t],nrow=M))>0,l]=1
      }
    }
    z=known.matrix
    # r=array(0,dim=dim(z))
    #turn on z's with caps on either side so we know they were in pop
    for(i in 1:n){
      idx=which(z[i,]==1)
      z[i,min(idx):max(idx)]=1
    }
    #turn on some augmented guys for z1
    # z[(n+1):M,1]=rbinom(M-n,1,psi)#add augmented guys to t=1 with psi.. not enough
    z1deal=psi*M-sum(z[,1])
    z[sample((n+1):M,z1deal),1]=1
    # r[,1]=z[,1]
    a=matrix(1,nrow=M,ncol=t) #a is available to be recruited
    # a[,1]=1-z[,1]
    a[which(z[,1]==1),]=0
    if(length(gamma)==(t-1)){
      gammasim=gamma
    }else{
      gammasim=rep(gamma,t-1)
    }
    if(length(phi)==(t-1)){
      phisim=phi
    }else{
      phisim=rep(phi,t-1)
    }

    for(i in 2:t){
      #Recruitment
      nrecruit=round(sum(z[,i-1])*gammasim[i-1])
      recruits=which(z[1:n,i-1]==0&z[1:n,i]==1)#who is already recruited based on y constraints?
      a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleft=nrecruit-length(recruits)
      if(nleft>0){
        cands=setdiff(which(a[,i-1]==1),recruits)#Who is available to recruit?
        cands=cands[!cands%in%1:n]#Don't mess with known guys
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        pick=sample(cands,nleft)
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Should probably raise M")
      }
      #Survival
      nlive=rbinom(1,sum(z[,i-1]),phisim[i-1])#How many to live
      livers=which(z[1:n,i-1]==1&z[1:n,i]==1)
      nleft=nlive-length(livers)
      a[livers,i]=0
      if(nleft>0){
        cands=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1)
        cands=cands[!cands%in%1:n]
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        pick=sample(cands,nleft)
        z[pick,i]=1
        a[pick,i:t]=0
      }
    }
    N=colSums(z)

    ll.z=matrix(0, M, t)
    Ez=Ez.cand=matrix(NA, M, t-1)
    ll.z[,1]=dbinom(z[,1], 1, psi, log=TRUE)
    gamma.prime <- numeric(t-1)
    if(length(gamma)==1){
      gammause=rep(gamma,t-1)
    }else{
      gammause=gamma
    }
    if(length(phi)==1){
      phiuse=rep(phi,t-1)
    }else{
      phiuse=phi
    }
    for(l in 2:t) {
      gamma.prime[l-1]=N[l-1]*gammause[l-1] / sum(a[,l-1])
      Ez[,l-1]=z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
      ll.z[,l]=dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
    }
    ll.z.cand=ll.z
    #Optimize starting locations given where they are trapped. Initalizing s1 and s2 at same locs
    s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(known.vector==1) #switch for those actually caught
    for(i in idx){
      trps=matrix(0,nrow=0,ncol=2)
      for(j in 1:t){ #loop over t to get all cap locs
        trps<- rbind(trps,X[[j]][y[i,,j]>0,1:2])
      }
      s1[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }

    #check to make sure everyone is in polygon
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
    }else{
      useverts=FALSE
    }
    if(useverts==TRUE){
      inside=rep(NA,nrow(s1))
      for(i in 1:nrow(s)){
        inside[i]=inout(s1[i,],vertices)
      }
      idx=which(inside==FALSE)
      if(length(idx)>0){
        for(i in 1:length(idx)){
          while(inside[idx[i]]==FALSE){
            s1[idx[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside[idx[i]]=inout(s1[idx[i],],vertices)
          }
        }
      }
    }
    #Initialize s2
    #Assign s2 to be s1 for all occasions
    s2=array(NA,dim=c(M,t,2))
    for(l in 1:t){
      s2[,l,]=s1
    }
    if(metamu){
      #update s2s for guys captured each year and add noise for uncaptured guys. More consistent with sigma_t>0
      for(l in 1:t){
        idx=which(rowSums(y[,,l])>0) #switch for those actually caught
        for(i in 1:M){
          if(i%in%idx){
            trps<- X[[l]][y[i,,l]>0,1:2]
            s2[i,l,]<- c(mean(trps[,1]),mean(trps[,2]))
          }else{
            inside=FALSE
            while(inside==FALSE){
              s2[i,l,]=c(rnorm(1,s1[i,1],sigma_t),rnorm(1,s1[i,2],sigma_t))
              inside=inout(s2[i,l,],vertices)
            }
          }
        }
      }
      ll.s2=dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE)
    }else{
      proppars$s1x=proppars$s1y=0
      ll.s2=matrix(NA,nrow=M,ncol=t)
      sigma_t=0
      proppars$sigma_t=0
    }
    if(niter<(nburn)){
      stop("niter is smaller than nburn")
    }
    D=lamd=ll.y=ll.y.cand=array(NA,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(i in 1:t){
      D[,1:nrow(X[[i]]),i]=e2dist(s2[,i,],X[[i]])
      if(length(lam0)==t&length(sigma)==t){
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam0)==1&length(sigma)==t){
        lamd[,,i]=lam0*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam0)==t&length(sigma)==1){
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma*sigma))
      }else{
        lamd[,,i]=lam0*exp(-D[,,i]^2/(2*sigma*sigma))
      }
    }
    #Calculate ll for observation model so we can check for -Inf values
    pd=pd.cand=1-exp(-lamd)
    for(l in 1:t){
      ll.y[,,l]= dbinom(y[,,l],K[l],pd[,,l]*z[,l],log=TRUE)
    }
    ll.y.cand=ll.y
    ll.y.t.sum=ll.y.cand.t.sum=apply(ll.y,3,sum) #ll summed for each year
    ll.y.sum=sum(ll.y.t.sum) #full ll sum
    #Check ll.y.sum for inf
    if(!is.finite(ll.y.sum)){
      stop("Detection function starting values produce -Inf log likelihood values. Try increasing sigma and/or lam0")
    }
    Xidx=unlist(lapply(X,dim))[seq(1,2*length(X)-1,2)]
    Xcpp=array(NA,dim=c(t,max(Xidx),2))
    for(l in 1:t){
      for(j in 1:Xidx[l]){
        Xcpp[l,j,]=as.numeric(X[[l]][j,])
      }
    }
    each=unlist(lapply(inits,length))[1:4]
    npar=sum(each)+t
    if(metamu){
      npar=npar+1
    }
    store=mcmc_Open(  lam0,  sigma,  gamma, gamma.prime,  phi, D,lamd, y, z, a,  s1, s2,
                      metamu, useverts, vertices, xlim, ylim, known.matrix, Xidx, Xcpp, K, Ez,  psi,
                      N, proppars$lam0, proppars$sigma, proppars$propz,  proppars$gamma, proppars$s1x,  proppars$s1y,
                      proppars$s2x,proppars$s2y,proppars$sigma_t,sigma_t,niter,nburn,nthin,npar,each)

    out=store[[1]]
    s1xout=store[[2]]
    s1yout=store[[3]]
    s2xout=store[[4]]
    s2yout=store[[5]]
    zout=store[[6]]
    # name out
    if(length(lam0)==t){
      lam0names=paste("lam0",1:t,sep="")
    }else{
      lam0names="lam0"
    }
    if(length(sigma)==t){
      sigmanames=paste("sigma",1:t,sep="")
    }else{
      sigmanames="sigma"
    }
    if(length(gamma)==(t-1)){
      gammanames=paste("gamma",1:(t-1),sep="")
    }else{
      gammanames="gamma"
    }
    if(length(phi)==(t-1)){
      phinames=paste("phi",1:(t-1),sep="")
    }else{
      phinames="phi"
    }
    Nnames=paste("N",1:t,sep="")
    if(metamu){
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"sigma_t")
    }else{
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames)
    }
    if(keepACs==TRUE){
      if(metamu){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout)
      }
    }else{
      list(out=out)
    }
  }
