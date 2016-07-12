mcmc.2side.ind2 <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam01=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    ###
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

    #Bernoulli Likelihood function
    func<- function(lamdL,lamdR,y.left,y.right,K,zL,zR,X){
      #convert lamd to pd (gaussian hazard model)
      pdL=1-exp(-lamdL)
      pdR=1-exp(-lamdR)
      #If data is M x K or 2 x K
      if(is.matrix(y.left)){
        trapno=matrix(rep(X[,3],M),nrow=M,byrow=TRUE) #trap number multiplier for left and right captures
        ones=trapno==1
        twos=trapno==2
        v1 <-  dbinom(y.left,K,ones*pdL+twos*(2*pdL-pdL*pdL),log=TRUE)
        v2 <-  dbinom(y.right,K,ones*pdR+twos*(2*pdR-pdR*pdR),log=TRUE)
        v1[zL==0,]<- 0
        v2[zR==0,]<- 0
        v=v1+v2
      }else{
        #If data is 1 x K
        ones=X[,3]==1
        twos=X[,3]==2
        v1 <-  zL*dbinom(y.left,K,ones*pdL+twos*(2*pdL-pdL*pdL),log=TRUE)
        v2 <-  zR*dbinom(y.right,K,ones*pdR+twos*(2*pdR-pdR*pdR),log=TRUE)
        v=v1+v2
      }
      v
    }

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam01","sigma","N1","N2"))
    sLxout<- sLyout<-sRxout<- sRyout<- zLout<-zRout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    DL<- e2dist(sL, X)
    DR<- e2dist(sR, X)
    lamdL<- lam01*exp(-DL*DL/(2*sigma*sigma))
    lamdR<- lam01*exp(-DR*DR/(2*sigma*sigma))

    for(i in 1:niter){
      #Update lam01
      lik.curr<-  sum( func(lamdL,lamdR, y.left,y.right,K,zL,zR,X) )
      lam01.cand<- rnorm(1,lam01,proppars$lam01)
      if(lam01.cand > 0){
        lamdLcand<- lam01.cand*exp(-DL*DL/(2*sigma*sigma))
        lamdRcand<- lam01.cand*exp(-DR*DR/(2*sigma*sigma))
        lik.new<-  sum( func(lamdLcand,lamdRcand,y.left,y.right,K,zL,zR,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          lam01<- lam01.cand
          lamdL=lamdLcand
          lamdR=lamdRcand
          lik.curr<- lik.new
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamdLcand<- lam01*exp(-DL*DL/(2*sigma.cand*sigma.cand))
        lamdRcand<- lam01*exp(-DR*DR/(2*sigma.cand*sigma.cand))
        lik.new<-  sum( func(lamdLcand,lamdRcand,y.left,y.right,K,zL,zR,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          sigma<- sigma.cand
          lamdL=lamdLcand
          lamdR=lamdRcand
          lik.curr<- lik.new
        }
      }
      #Update psi gibbs
      ## probability of not being captured in a trap AT ALL
      trapno=matrix(rep(X[,3],M),nrow=M,byrow=TRUE) #trap number multiplier for left and right captures
      ones=trapno==1
      twos=trapno==2
      #left data set
      pdL=1-exp(-lamdL)
      pdLtraps=ones*pdL+twos*(2*pdL-pdL*pdL)
      pbarL=(1-pdLtraps)^K
      prob0L<- exp(rowSums(log(pbarL)))
      fcL<- prob0L*psi1/(prob0L*psi1 + 1-psi1)
      zL[cap.vector1==0]<- rbinom(sum(cap.vector1 ==0), 1, fcL[cap.vector1==0])
      #left data set
      pdR=1-exp(-lamdR)
      pdRtraps=ones*pdR+twos*(2*pdR-pdR*pdR)
      pbarR=(1-pdRtraps)^K
      prob0R<- exp(rowSums(log(pbarR)))
      fcR<- prob0R*psi2/(prob0R*psi2 + 1-psi2)
      zR[cap.vector2==0]<- rbinom(sum(cap.vector2 ==0), 1, fcR[cap.vector2==0])
      lik.curr<-  sum( func(lamdL,lamdR, y.left,y.right,K,zL,zR,X) )
      psi1 <- rbeta(1, 1 + sum(zL), 1 + M - sum(zL))
      psi2 <- rbeta(1, 1 + sum(zR), 1 + M - sum(zR))

      # Now we have to update the activity centers
      #left guys first
      for (j in 1:M) {
        Scand <- c(rnorm(1, sL[j, 1], proppars$sx), rnorm(1, sL[j, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamdLthisj<- lam01*exp(-DL[j,]*DL[j,]/(2*sigma*sigma))
          lamdRthisj<- lam01*exp(-DR[j,]*DR[j,]/(2*sigma*sigma))
          lamdLcand<-lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamdLthisj,lamdRthisj,y.left[j,],y.right[j,],K,zL[j],zR[j],X))
          llcand<- sum(func(lamdLcand,lamdRthisj,y.left[j,],y.right[j,],K,zL[j],zR[j],X))
          if (runif(1) < exp(llcand - llS)) {
            sL[j, ] <- Scand
            DL[j, ] <- dtmp
            lamdL[j, ] <- lamdLcand
          }
        }
      }
      #right guys next
      for (j in 1:M) {
        Scand <- c(rnorm(1, sR[j, 1], .2), rnorm(1, sR[j, 2], .2))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamdLthisj<- lam01*exp(-DL[j,]*DL[j,]/(2*sigma*sigma))
          lamdRthisj<- lam01*exp(-DR[j,]*DR[j,]/(2*sigma*sigma))
          lamdRcand<-lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamdLthisj,lamdRthisj,y.left[j,],y.right[j,],K,zL[j],zR[j],X))
          llcand<- sum(func(lamdLthisj,lamdRcand,y.left[j,],y.right[j,],K,zL[j],zR[j],X))
          if (runif(1) < exp(llcand - llS)) {
            sR[j, ] <- Scand
            DR[j, ] <- dtmp
            lamdR[j, ] <- lamdRcand
          }
        }
      }
      #Do we record output on this iteration?
      if(i>nburn&i%%nthin==0){
        sLxout[idx,]<- sL[,1]
        sLyout[idx,]<- sL[,2]
        sRxout[idx,]<- sR[,1]
        sRyout[idx,]<- sR[,2]
        zLout[idx,]<- zL
        zRout[idx,]<- zR
        out[idx,]<- c(lam01,sigma ,sum(zL),sum(zR))
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      list(out=out, sLxout=sLxout, sLyout=sLyout,sRxout=sRxout, sRyout=sRyout, zLout=zLout,zRout=zRout)
    }else{
      list(out=out)
    }
  }
