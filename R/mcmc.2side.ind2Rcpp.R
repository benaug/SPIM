mcmc.2side.ind2Rcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam01=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=FALSE){
    library(abind)
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
    lam01<- inits$lam01
    sigma<- inits$sigma

    #augment both data sets
    left<- abind(left,array(0, dim=c( M-dim(left)[1],3, K, J)), along=1)
    right<-abind(right,array(0,dim=c( M-dim(right)[1],3,K,J)), along=1)
    y.left<- apply(left[,2,,], c(1,3), sum)
    y.right<- apply(right[,3,,], c(1,3), sum)

    #Optimize starting locations given where they are trapped.
    sL<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    sR<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    if(useverts==TRUE){
      inside1=rep(NA,nrow(sL))
      inside2=rep(NA,nrow(sR))
      for(i in 1:nrow(sL)){
        inside1[i]=inout(sL[i,],vertices)
      }
      for(i in 1:nrow(sR)){
        inside2[i]=inout(sR[i,],vertices)
      }
      idx1=which(inside1==FALSE)
      if(length(idx1)>0){
        for(i in 1:length(idx1)){
          while(inside[idx1[i]]==FALSE){
            sL[idx1[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside1[idx1[i]]=inout(sL[idx1[i],],vertices)
          }
        }
      }
      idx2=which(inside1==FALSE)
      if(length(idx2)>0){
        for(i in 1:length(idx2)){
          while(inside[idx2[i]]==FALSE){
            sR[idx2[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside2[idx2[i]]=inout(sR[idx2[i],],vertices)
          }
        }
      }
    }
    cap.vector1<- (rowSums(y.left)>0)*1
    cap.vector2<- (rowSums(y.right)>0)*1
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

    zL=zR=rep(0,M)
    zL[cap.vector1==1]=1
    zR[cap.vector2==1]=1
    zL[sample(which(zL==0),sum(zL==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
    zR[sample(which(zR==0),sum(zR==0)/2)]=1

    DL<- e2dist(sL, X)
    DR<- e2dist(sR, X)
    lamdL<- lam01*exp(-DL*DL/(2*sigma*sigma))
    lamdR<- lam01*exp(-DR*DR/(2*sigma*sigma))


    #Run MCMC
    store=MCMCInd2(lam01,sigma,y.left,y.right,zL,zR,X,K,DL,DR,cap.vector1,cap.vector2,sL,sR,psi1,psi2,xlim,ylim,useverts,vertices,
                   proppars$lam01,proppars$sigma,proppars$sx,proppars$sy,niter,nburn,nthin)

    if(keepACs){
      list(out=store[[1]], sxout=store[[2]], syout=store[[3]], zout=store[[4]],data=data)
    }else{
      list(out=store[[1]],data=data)
    }
  }
