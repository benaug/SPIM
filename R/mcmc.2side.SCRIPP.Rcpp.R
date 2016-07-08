#' Spatial partial identity MCMC algorithm with inhomogenous density using Rcpp.
#' @param data a list produced by sim2side or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam01=lam01,lam02=lam02,sigma=sigma)
#' @param  swap number of IDs to swap on each MCMC iteration
#' @param  swap.tol the search radius within which to search for partial ID activity centers to match with
#' @param proppars a list of tuning parameters for the proposal distributions
#' @return  a list with the posteriors for the SCR parameters (out), s, z, ID_L and ID_R
#' @author Ben Augustine, Andy Royle
#' @description This function runs the MCMC algorithm for the spatial partial identity model.  The data list should have the following elements:
#' 1.  both, a n_both x 3 x K x J both side data array.  If n_both=0 as in an all single camera study, the first dimension is 0 and the data
#' set should still have 4 dimensions, 0 x 3 x K x J.
#' 2.  left, a n_left x 3 x K x J left side data array.
#' 3.  right, a n_right x 3 x K x J left side data array.
#' 4.  IDknown a vector listing the index of complete identity individuals.  It is assumed individuals are sorted such that the complete identity
#' individuals are listed first.  So if there are 7 complete identity individuals, IDknown=1:7.
#' 5. X a matrix with the X and Y trap locations in the first two columns and the number of cameras (1 or 2) at each trap in the third.
#' 6. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' @examples
#' \dontrun{N=50
#'p01=0.13
#'p02=0.2
#'lam01=-log(1-p01)
#'lam02=-log(1-p02)
#'sigma=0.50
#'K=5
#'buff=2
#'niter=1000 #should run more than this and discard a burn in
#'nburn=1
#'nthin=1
#'xlim<- c(1,10)
#'ylim<- c(1,10)
#'X<- expand.grid(3:8,3:8) #6x6 trapping array
#'X=cbind(X,1) #add number of cameras at each trap
#'X[which(X[,2]%in%c(4,7)),3]=2 #switch the second and fifth row of traps to double cameras
#'#Simulate some data
#'data=sim2side(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff)
#'inits=list(psi=0.5,lam01=lam01,lam02=lam02,sigma=sigma)
#'a=Sys.time()
#'store=mcmc.2side(data,niter=niter,nburn=nburn,nthin=nthin, M = 100,inits=inits,swap=10)
#'b=Sys.time()
#'b-a
#'#plot posteriors
#'par(mfrow=c(2,2))
#'plot(store$out[,1],main="lam01",type="l")
#'plot(store$out[,2],main="lam02",type="l")
#'plot(store$out[,3],main="sigma",type="l")
#'plot(store$out[,4],main="N",type="l")
#'par(mfrow=c(1,1))}
#' @export
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


    store=SCRRcpp::MCMC2( lam01, lam02,  sigma,beta0,beta1,y.both,  y.left,  y.right,  z,  X, K,D, Nfixed,known.vector,
                          ID_L,  ID_R, swap, swap.tol,aperm(left[,2,,],c(2,3,1)),aperm(right[,3,,],c(2,3,1)),
                          as.matrix(s[,2:3]),s[,1],psi,grid,cellArea,EN,proppars$lam01,proppars$lam02,proppars$sigma,proppars$beta0,
                          proppars$beta1,proppars$sx,proppars$sy,niter,nburn,nthin,updates)
    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]], ID_L=store[[4]],ID_R=store[[5]],zout=store[[6]])
    }else{
      list(out=store[[1]])
    }
  }
