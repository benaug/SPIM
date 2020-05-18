mcmc.2sidetfFullRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,swap=10,swap.tol=1,
           proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),storeLatent=FALSE){
    library(abind)
    y.both<-data$both
    y.left.obs<-data$left
    y.right.obs<-data$right
    if(length(dim(y.both))==3){
      y.both=apply(y.both,c(1,2),sum)
    }
    if(length(dim(y.left.obs))==2){
      stop("Must use 3D single side data with 2D trap file")
    }
    if(dim(y.right.obs)[1]>dim(y.left.obs)[1]){
      storeL=y.left.obs
      storeR=y.right.obs
      dimL=dim(y.left.obs)
      y.left.obs=array(0,dim=dim(y.right.obs))
      y.right.obs=array(0,dim=dimL)
      y.left.obs=storeR
      y.right.obs=storeL
      warning("Right side data set larger than left so I switched them for convenience")
    }
    if(!("tf"%in%names(data))){
      stop("This function requires a 2D trap file. Tell Ben there is an error")
    }
    if(!is.matrix(data$tf)){
      stop("This is the wrong function for a 2D trap file. Tell Ben something went wrong")
    }
    tf1=data$tf
    tf2=rowSums(tf1==2)
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- data$K
    IDknown<- data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(y.left.obs)[1]-Nfixed
    nright<-dim(y.right.obs)[1]-Nfixed
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
      xlim=c(min(vertices[,1]),max(vertices[,1]))
      ylim=c(min(vertices[,2]),max(vertices[,2]))
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
    sigma<- inits$sigma
    lam01<- inits$lam01
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
    
    #augment 3 data sets
    y.both<- abind(y.both,array(0, dim=c( M-dim(y.both)[1],J)), along=1)
    y.left.obs<- abind(y.left.obs,array(0, dim=c( M-dim(y.left.obs)[1],J,K)), along=1)
    y.right.obs<-abind(y.right.obs,array(0,dim=c( M-dim(y.right.obs)[1],J,K)), along=1)
    
    #sort to minimize distance between initial matches. Skip if no single sides.
    # if(nleft>0|nright>0){
    if(nleft>0&nright>0){ #change 1/14/19
      y.left.tmp=apply(y.left.obs,c(1,2),sum)
      y.right.tmp=apply(y.right.obs,c(1,2),sum)
      IDs<- LRmatch(M=M,left=y.left.tmp, nleft=nleft, right=y.right.tmp, nright=nright, X, Nfixed=Nfixed)
      #Add unused augmented indivuals back in
      notusedL<- (1:M)[is.na(match(1:M,IDs$ID_L))]
      ID_L<-c(IDs$ID_L,notusedL)
      notusedR<- (1:M)[is.na(match(1:M,IDs$ID_R))]
      ID_R<-c(IDs$ID_R,notusedR)
    }else{
      ID_R=ID_L=1:M
    }
    
    #reoder left and right possibly true data sets
    y.left.true<- y.left.obs[order(ID_L),,]
    y.right.true<- y.right.obs[order(ID_R),,]
    
    #Make initial complete data set
    tmpdata<- y.both + apply(y.left.true+y.right.true,c(1,2),sum)
    z=1*(apply(tmpdata,1,sum)>0)
    z[sample(which(z==0),sum(z==0)*psi)]=1 #switch some uncaptured z's to 1.
    
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
    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )
    D<- e2dist(s, X)
    lamd1<- lam01*exp(-D*D/(2*sigma*sigma))
    lamd2<- lam02*exp(-D*D/(2*sigma*sigma))
    #check if ll finite
    tf1.tmp=abind(tf1,tf1,along=0)
    for(i in 3:M){
      tf1.tmp=abind(tf1.tmp,tf1,along=1)
    }
    tf2.tmp=matrix(rep(tf2,M),nrow=M,byrow=TRUE)
    
    pd1=1-exp(-lamd1)
    pd12=2*pd1-pd1*pd1
    pd2=1-exp(-lamd2)
    lamd1.cand=lamd1
    lamd2.cand=lamd2
    pd1.cand=pd1
    pd2.cand=pd2
    zero.guys<- apply(y.both,1,sum)+apply(y.left.true + y.right.true ,1,sum) == 0
    ones=1*(tf1.tmp==1)
    twos=1*(tf1.tmp==2)
    sumones=apply(ones,c(1,2),sum)
    sumtwos=apply(twos,c(1,2),sum)
    ll.y.both=dbinom(y.both,tf2.tmp,z*pd2,log=TRUE)
    ll.y.left=dbinom(y.left.true*ones,ones,z*pd1,log=TRUE)+
      dbinom(y.left.true*twos,twos,z*pd12,log=TRUE)
    ll.y.right=dbinom(y.right.true*ones,ones,z*pd1,log=TRUE)+
      dbinom(y.right.true*twos,twos,z*pd12,log=TRUE)
    if(!is.finite(sum(ll.y.both)))stop("Both side likelihood not finite. Make sure all camera stations recording both side captures have 2 cameras. Then try changing lam02 or sigma inits.")
    if(!is.finite(sum(ll.y.left)))stop("Left side likelihood not finite. Try changing lam01 or sigma inits.")
    if(!is.finite(sum(ll.y.right)))stop("right side likelihood not finite. Try changing lam01 or sigma inits.")
    
    #Run MCMC
    store=MCMCtf2(lam01,lam02,sigma,lamd1,lamd2,y.both, y.left.true,y.right.true,
                  y.left.obs,y.right.obs,
                  z,X, tf1,tf2,D,Nfixed,known.vector,ID_L,ID_R,swap,swap.tol,
                  s,psi,xlim,ylim,useverts,vertices,proppars$lam01,proppars$lam02,proppars$sigma,proppars$sx,
                  proppars$sy,niter,nburn,nthin,updates,storeLatent=storeLatent)
    out=store[[1]]
    colnames(out)=c("lam01","lam02","sigma","N","psi")
    if(storeLatent){
      list(out=out, sxout=store[[2]], syout=store[[3]], ID_Lout=store[[4]],ID_Rout=store[[5]],zout=store[[6]])
    }else{
      list(out=out)
    }
  }
