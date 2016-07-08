#' Run MCMC algorithm for "heuristic estimator" with both left and right data set
#' @param data a list produced by sim2side or in the same format
#' @param niter
#' @param  nburn
#' @param nthin
#' @param M
#' @param inits
#' @return  a list with output
#' @author Ben Augustine, Andy Royle
#' @description This function runs the MCMC algorithm for the naive independence estimator when there are 3 data sets (both, left, and right).
#' The data list should have the following elements:
#' 1.  both, a n_both x 3 x K x J both side data array.  If n_both=0 as in an all single camera study, the first dimension is 0 and the data
#' set should still have 4 dimensions, 0 x 3 x K x J.
#' 2.  left, a n_left x 3 x K x J left side data array.
#' 3.  right, a n_right x 3 x K x J left side data array.
#' 4. X a matrix with the X and Y trap locations in the first two columns and the number of cameras (1 or 2) at each trap in the third.
#' 5. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' @export

mcmc.2side.ind3 <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2)){
    ###
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
    zB[sample(which(zB==0),sum(zB==0)/2)]=1
    zL[sample(which(zL==0),sum(zL==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
    zR[sample(which(zR==0),sum(zR==0)/2)]=1

    #Bernoulli Likelihood function
    func<- function(lamdB,lamdL,lamdR,y.both,y.left,y.right,K,zB,zL,zR,X){
      #convert lamd to pd (gaussian hazard model)
      pdB=1-exp(-lamdB)
      pdL=1-exp(-lamdL)
      pdR=1-exp(-lamdR)
      #If data is M x K or 2 x K
      if(is.matrix(y.left)){
        v0 <- dbinom(y.both,K,pdB,log=TRUE)
        v0[,X[,3]==1]=0 #cancel out contributions from single traps
        trapno=matrix(rep(X[,3],M),nrow=M,byrow=TRUE) #trap number multiplier for left and right captures
        ones=trapno==1
        twos=trapno==2
        v1 <-  dbinom(y.left,K,ones*pdL+twos*(2*pdL-pdL*pdL),log=TRUE)
        v2 <-  dbinom(y.right,K,ones*pdR+twos*(2*pdR-pdR*pdR),log=TRUE)
        v0[zB==0,]=0
        v1[zL==0,]<- 0
        v2[zR==0,]<- 0
        v=v0+v1+v2
      }else{
        #If data is 1 x K
        v0 <- zB*dbinom(y.both,K,pdB,log=TRUE)
        v0[X[,3]==1]=0 #cancel out contributions from single traps
        ones=X[,3]==1
        twos=X[,3]==2
        v1 <-  zL*dbinom(y.left,K,ones*pdL+twos*(2*pdL-pdL*pdL),log=TRUE)
        v2 <-  zR*dbinom(y.right,K,ones*pdR+twos*(2*pdR-pdR*pdR),log=TRUE)
        v=v0+v1+v2
      }
      v
    }

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

    for(i in 1:niter){
      #Update lam01
      lik.curr<-  sum( func(lamdB,lamdL,lamdR, y.both,y.left,y.right,K,zB,zL,zR,X) )
      lam01.cand<- rnorm(1,lam01,proppars$lam01)
      if(lam01.cand > 0){
        lamdLcand<- lam01.cand*exp(-DL*DL/(2*sigma*sigma))
        lamdRcand<- lam01.cand*exp(-DR*DR/(2*sigma*sigma))
        lik.new<-  sum( func(lamdB,lamdLcand,lamdRcand,y.both,y.left,y.right,K,zB,zL,zR,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          lam01<- lam01.cand
          lamdL=lamdLcand
          lamdR=lamdRcand
          lik.curr<- lik.new
        }
      }#Update lam02
      lam02.cand<- rnorm(1,lam02,proppars$lam02)
      if(lam02.cand > 0){
        lamdBcand<- lam02.cand*exp(-DB*DB/(2*sigma*sigma))
        lik.new<-  sum( func(lamdBcand,lamdL,lamdR,y.both,y.left,y.right,K,zB,zL,zR,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          lam02<- lam02.cand
          lamdB=lamdBcand
          lik.curr<- lik.new
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamdBcand<- lam02*exp(-DB*DB/(2*sigma.cand*sigma.cand))
        lamdLcand<- lam01*exp(-DL*DL/(2*sigma.cand*sigma.cand))
        lamdRcand<- lam01*exp(-DR*DR/(2*sigma.cand*sigma.cand))
        lik.new<-  sum( func(lamdBcand,lamdLcand,lamdRcand,y.both,y.left,y.right,K,zB,zL,zR,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          sigma<- sigma.cand
          lamdB=lamdBcand
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
      #both data set
      pdB=1-exp(-lamdB)
      pbarB=(1-pdB)^K
      prob0B<- exp(rowSums(log(pbarB)))
      fcB<- prob0B*psi1/(prob0B*psi1 + 1-psi1)
      zB[cap.vector1==0]<- rbinom(sum(cap.vector1 ==0), 1, fcB[cap.vector1==0])
      #left data set
      pdL=1-exp(-lamdL)
      pdLtraps=ones*pdL+twos*(2*pdL-pdL*pdL)
      pbarL=(1-pdLtraps)^K
      prob0L<- exp(rowSums(log(pbarL)))
      fcL<- prob0L*psi2/(prob0L*psi2 + 1-psi2)
      zL[cap.vector2==0]<- rbinom(sum(cap.vector2 ==0), 1, fcL[cap.vector2==0])
      #right data set
      pdR=1-exp(-lamdR)
      pdRtraps=ones*pdR+twos*(2*pdR-pdR*pdR)
      pbarR=(1-pdRtraps)^K
      prob0R<- exp(rowSums(log(pbarR)))
      fcR<- prob0R*psi3/(prob0R*psi3 + 1-psi3)
      zR[cap.vector3==0]<- rbinom(sum(cap.vector3==0), 1, fcR[cap.vector3==0])
      lik.curr<-  sum( func(lamdB,lamdL,lamdR,y.both,y.left,y.right,K,zB,zL,zR,X) )
      psi1 <- rbeta(1, 1 + sum(zB), 1 + M - sum(zB))
      psi2 <- rbeta(1, 1 + sum(zL), 1 + M - sum(zL))
      psi3 <- rbeta(1, 1 + sum(zR), 1 + M - sum(zR))

      # Now we have to update the activity centers
      #both guys first
      for (j in 1:M) {
        Scand <- c(rnorm(1, sB[j, 1], proppars$sx), rnorm(1, sB[j, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamdBthisj<- lam02*exp(-DB[j,]*DB[j,]/(2*sigma*sigma))
          lamdLthisj<- lam01*exp(-DL[j,]*DL[j,]/(2*sigma*sigma))
          lamdRthisj<- lam01*exp(-DR[j,]*DR[j,]/(2*sigma*sigma))
          lamdBcand<-lam02*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamdBthisj,lamdLthisj,lamdRthisj,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
          llcand<- sum(func(lamdBcand,lamdLthisj,lamdRthisj,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
          if (runif(1) < exp(llcand - llS)) {
            sB[j, ] <- Scand
            DB[j, ] <- dtmp
            lamdB[j, ] <- lamdBcand
          }
        }
      }
      #left guys next
      for (j in 1:M) {
        Scand <- c(rnorm(1, sL[j, 1], proppars$sx), rnorm(1, sL[j, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamdBthisj<- lam02*exp(-DB[j,]*DB[j,]/(2*sigma*sigma))
          lamdLthisj<- lam01*exp(-DL[j,]*DL[j,]/(2*sigma*sigma))
          lamdRthisj<- lam01*exp(-DR[j,]*DR[j,]/(2*sigma*sigma))
          lamdLcand<-lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamdBthisj,lamdLthisj,lamdRthisj,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
          llcand<- sum(func(lamdBthisj,lamdLcand,lamdRthisj,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
          if (runif(1) < exp(llcand - llS)) {
            sL[j, ] <- Scand
            DL[j, ] <- dtmp
            lamdL[j, ] <- lamdLcand
          }
        }
      }
      #right guys next
      for (j in 1:M) {
        Scand <- c(rnorm(1, sR[j, 1], proppars$sx), rnorm(1, sR[j, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamdBthisj<- lam02*exp(-DB[j,]*DB[j,]/(2*sigma*sigma))
          lamdLthisj<- lam01*exp(-DL[j,]*DL[j,]/(2*sigma*sigma))
          lamdRthisj<- lam01*exp(-DR[j,]*DR[j,]/(2*sigma*sigma))
          lamdRcand<-lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamdBthisj,lamdLthisj,lamdRthisj,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
          llcand<- sum(func(lamdBthisj,lamdLthisj,lamdRcand,y.both[j,],y.left[j,],y.right[j,],K,zB[j],zL[j],zR[j],X))
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
        out[idx,]<- c(lam01,lam02,sigma ,(sum(zB)+sum(zL)+sum(zR))/3)
        idx=idx+1
      }
    }  # end of MCMC algorithm

    list(out=out, sLxout=sLxout, sLyout=sLyout,sRxout=sRxout, sRyout=sRyout, zLout=zLout,zRout=zRout)
  }
