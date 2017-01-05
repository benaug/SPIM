mcmc.2sideRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,swap=10,swap.tol=1,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=FALSE){
    library(abind)
    both<-data$both
    left<-data$left
    right<-data$right
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
    #If using polygon state space
    nright<-dim(right)[1]-Nfixed
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
      ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
      vertices=rbind(xlim,ylim)
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    ##pull out initial values
    psi<- inits$psi
    lam01<- inits$lam01
    lam02<- inits$lam02
    sigma<- inits$sigma

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

    #sort to minimize distance between initial matches. Skip if no single sides. Need to fix if only lefts or only rights >Nfixed
    if(nleft>0&nright>0){
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
    z=1*(apply(tmpdata,1,sum)>0)
    z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?

    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(tmpdata)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[tmpdata[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }
    if(useverts==TRUE){
      inside=rep(NA,nrow(s))
      for(i in 1:nrow(s)){
        inside[i]=SCRRcpp::inout(s[i,],vertices)
      }
      idx=which(inside==FALSE)
      if(length(idx)>0){
        for(i in 1:length(idx)){
          while(inside[idx[i]]==FALSE){
            s[idx[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside[idx[i]]=SCRRcpp::inout(s[idx[i],],vertices)
          }
        }
      }
    }
    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )

    #Create initial data sets
    y.both<- apply(both[,1,,], c(1,3), sum)
    y.left<- apply(left[order(ID_L),2,,], c(1,3), sum)
    y.right<- apply(right[order(ID_R),3,,], c(1,3), sum)
    D<- e2dist(s, X)
    lamd1<- lam01*exp(-D*D/(2*sigma*sigma))
    lamd2<- lam02*exp(-D*D/(2*sigma*sigma))
    zero.guys<- apply(y.both+y.left + y.right ,1,sum) == 0

    #Run MCMC
    store=SPIM::MCMC( lam01, lam02,  sigma,y.both,  y.left,  y.right,  z,  X, K,D, Nfixed,known.vector,
                         ID_L,  ID_R, swap, swap.tol,aperm(left[,2,,],c(2,3,1)),aperm(right[,3,,],c(2,3,1)),
                         s,psi,xlim,ylim,useverts,vertices,proppars$lam01,proppars$lam02,proppars$sigma,proppars$sx,
                         proppars$sy,niter,nburn,nthin,updates)
    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]], ID_Lout=store[[4]],ID_Rout=store[[5]],zout=store[[6]])
    }else{
      list(out=store[[1]])
    }
  }
