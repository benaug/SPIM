mcmc.2side.SCRIPP.Rcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, cellArea=cellArea,inits=inits,swap=10,swap.tol=1,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    library(abind)
    both<-data$both
    left<-data$left
    right<-data$right
    if(dim(right)[1]>dim(left)[1]){
      storeL=left[,2,,]
      storeR=right[,3,,]
      dimL=dim(left)
      left=array(0,dim=dim(right))
      right=array(0,dim=dimL)
      left[,2,,]=storeR
      right[,3,,]=storeL
      warning("Right side data set larger than left so I switched them for convenience")
    }
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(left)[3]
    IDknown<- data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(left)[1]-Nfixed
    nright<-dim(right)[1]-Nfixed
    ##pull out initial values
    psi<- inits$psi
    lam01<- inits$lam01
    sigma<- inits$sigma
    beta0<- inits$beta0
    beta1<- inits$beta1
    grid=as.matrix(data$grid)
    npix=nrow(grid)
    lam02<- inits$lam02
    #Figure out what needs to be updated
    uplam01=uplam02=upLeft=upRight=TRUE
    if(lam01==0){
      uplam01=FALSE
      upIDs=FALSE
    }
    if(all(X[,3]==1)|lam02==0){
      uplam02=FALSE
    }
    if(upLeft==TRUE&(nleft==0)){
      upLeft=FALSE
    }
    if(upRight==TRUE&(nright==0)){
      upRight=FALSE
    }
    updates=c(uplam01,uplam02,upLeft,upRight)

    #augment both data
    both<- abind(both,array(0, dim=c( M-dim(both)[1],3,K, J)), along=1)
    left<- abind(left,array(0, dim=c( M-dim(left)[1],3, K, J)), along=1)
    right<-abind(right,array(0,dim=c( M-dim(right)[1],3,K,J)), along=1)

    #sort to minimize distance between initial matches. Skip if no single sides.
    if(nleft>0|nright>0){
      IDs<- LRmatch(M=M,left=left, nleft=nleft, right=right, nright=nright, X, Nfixed=Nfixed)
      #Add unused augmented indivuals back in
      notusedL<- (1:M)[is.na(match(1:M,IDs$ID_L))]
      ID_L<-c(IDs$ID_L,notusedL)
      notusedR<- (1:M)[is.na(match(1:M,IDs$ID_R))]
      ID_R<-c(IDs$ID_R,notusedR)
    }else{
      ID_R=ID_L=1:M
    }

    #Make initial complete data set
    tmpdata<- both + left[order(ID_L),,,] + right[order(ID_R),,,]
    tmpdata<- apply(tmpdata,c(1,4),sum)
    mu=exp(beta0+beta1*1*(grid[,3]==2))*cellArea
    EN=sum(mu)
    psi <- EN/M
    if (psi > 1){
      stop("Bad initial values for beta0 or beta1. Or M is too low")
    }
    z=1*(apply(tmpdata,1,sum)>0)
    N=sum(z)
    N2=round(psi*EN)

    if(N<N2){
      z[sample(which(z==0),N2-N)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
      N=sum(z)
    }

    #Optimize starting locations given where they are trapped.
    cell=sample(1:npix,M)#assign random locations
    s<- cbind(cell,grid[cell,1:2])
    idx=which(rowSums(tmpdata)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[tmpdata[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      loc<- c(mean(trps[,1]),mean(trps[,2]))
      d=sqrt((loc[1]-grid[,1])^2+(loc[2]-grid[,2])^2)
      cell=which(d==min(d))
      if(length(cell)>1){
        cell=cell[1]
      }
      s[i,]=c(cell,grid[cell,1:2])
    }

    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=7)
    dimnames(out)<-list(NULL,c("lam01","lam02","sigma","N","EN","beta0","beta1"))
    sxout<- syout<- zout<- ID_Lout<- ID_Rout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    y.both<- apply(both[,1,,], c(1,3), sum)
    y.left<- apply(left[order(ID_L),2,,], c(1,3), sum)
    y.right<- apply(right[order(ID_R),3,,], c(1,3), sum)
    D<- e2dist(s[,2:3], X)


    store=SPIM::MCMC2( lam01, lam02,  sigma,beta0,beta1,y.both,  y.left,  y.right,  z,  X, K,D, Nfixed,known.vector,
                          ID_L,  ID_R, swap, swap.tol,aperm(left[,2,,],c(2,3,1)),aperm(right[,3,,],c(2,3,1)),
                          as.matrix(s[,2:3]),s[,1],psi,grid,cellArea,EN,proppars$lam01,proppars$lam02,proppars$sigma,proppars$beta0,
                          proppars$beta1,proppars$sx,proppars$sy,niter,nburn,nthin,updates)
    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]], ID_L=store[[4]],ID_R=store[[5]],zout=store[[6]])
    }else{
      list(out=store[[1]])
    }
  }
