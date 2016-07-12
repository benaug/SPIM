mcmc.2side.ind3Rcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=FALSE){
    library(abind)
    both<-data$both
    left<-data$left
    right<-data$right
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(left)[3]
    IDknown<- data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(data$left)[1]-Nfixed
    nright<-dim(data$right)[1]-Nfixed
    buff<- data$buff
    xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
    ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
    #If using polygon state space
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
    }else{
      useverts=FALSE
      vertices=rbind(xlim,ylim)
    }

    ##pull out initial values
    psi1<- inits$psi1
    psi2<- inits$psi2
    psi3<- inits$psi3
    lam01<- inits$lam01
    lam02<- inits$lam02
    sigma<- inits$sigma

    #augment both data
    both<- abind(both,array(0, dim=c( M-dim(both)[1],3,K, J)), along=1)
    left<- abind(left,array(0, dim=c( M-dim(left)[1],3, K, J)), along=1)
    right<-abind(right,array(0,dim=c( M-dim(right)[1],3,K,J)), along=1)

    #move both data set to left and right components
    y.both<- apply(both[,1,,], c(1,3), sum)
    y.left<- apply(left[,2,,], c(1,3), sum)
    y.right<- apply(right[,3,,], c(1,3), sum)

    #Optimize starting locations given where they are trapped.
    sL<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    sR<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    sB<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    if(useverts==TRUE){
      inside1=rep(NA,nrow(sL))
      inside2=rep(NA,nrow(sR))
      inside3=rep(NA,nrow(sB))
      for(i in 1:nrow(sL)){
        inside1[i]=inout(sL[i,],vertices)
      }
      for(i in 1:nrow(sR)){
        inside2[i]=inout(sR[i,],vertices)
      }
      for(i in 1:nrow(sB)){
        inside3[i]=inout(sB[i,],vertices)
      }
      idx1=which(inside1==FALSE)
      if(length(idx1)>0){
        for(i in 1:length(idx1)){
          while(inside1[idx1[i]]==FALSE){
            sL[idx1[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside1[idx1[i]]=inout(sL[idx1[i],],vertices)
          }
        }
      }
      idx2=which(inside2==FALSE)
      if(length(idx2)>0){
        for(i in 1:length(idx2)){
          while(inside[idx2[i]]==FALSE){
            sR[idx2[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside2[idx2[i]]=inout(sR[idx2[i],],vertices)
          }
        }
      }
      idx3=which(inside3==FALSE)
      if(length(idx3)>0){
        for(i in 1:length(idx3)){
          while(inside[idx3[i]]==FALSE){
            sR[idx3[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside3[idx3[i]]=inout(sR[idx3[i],],vertices)
          }
        }
      }
    }
    cap.vector1<- (rowSums(y.both)>0)*1
    cap.vector2<- (rowSums(y.left)>0)*1
    cap.vector3<- (rowSums(y.right)>0)*1
    #for captured guys, pick smart starting location
    idx=which(rowSums(y.left)>0)
    for(i in idx){
      trps<- X[y.left[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      sL[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }
    idx=which(rowSums(y.right)>0)
    for(i in idx){
      trps<- X[y.right[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      sR[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }
    idx=which(rowSums(y.both)>0)
    for(i in idx){
      trps<- X[y.both[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      sB[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }

    zL=zR=zB=rep(0,M)
    zB[cap.vector1==1]=1
    zL[cap.vector2==1]=1
    zR[cap.vector3==1]=1
    zB[sample(which(zB==0),sum(zB==0)/1.5)]=1
    zL[sample(which(zL==0),sum(zL==0)/1.5)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
    zR[sample(which(zR==0),sum(zR==0)/1.5)]=1


    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam01","lam02","sigma","N"))
    sBxout<- sByout<-sLxout<- sLyout<-sRxout<- sRyout<- zLout<-zRout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    DB<- e2dist(sB, X)
    DL<- e2dist(sL, X)
    DR<- e2dist(sR, X)
    lamdB<- lam02*exp(-DB*DB/(2*sigma*sigma))
    lamdL<- lam01*exp(-DL*DL/(2*sigma*sigma))
    lamdR<- lam01*exp(-DR*DR/(2*sigma*sigma))


    #Run MCMC
    store=MCMCInd3(lam01,lam02,sigma,y.both,y.left,y.right,zB,zL,zR,X,K,DB,DL,DR,cap.vector1,cap.vector2,cap.vector3,
                   sB,sL,sR,psi1,psi2,psi3,xlim,ylim,useverts,vertices,proppars$lam01,proppars$lam02,proppars$sigma,proppars$sx,proppars$sy,niter,nburn,nthin)

    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]], zout=store[[4]],data=data)
    }else{
      list(out=store[[1]],data=data)
    }
  }
