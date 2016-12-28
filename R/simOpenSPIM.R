e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


cellprobsOpen<- function(lamd){
  # For gaussian hazard model convert lamda(s,x) to p(s,x)
  N<- dim(lamd)[1]
  J<- dim(lamd)[2] # traps
  t<- dim(lamd)[3]
  # just 2 outcomes: left captures and right captures
  pmat<- array(NA,dim=c(N,J,t,2))
  for(l in 1:t){
    for(j in 1:J){
      pmat[,j,l,]<- 1-exp(-lamd[,j,l,])
    }
  }
  pmat
}

#' Simulate data from open population SPIM.  Update help file later.  Add polygon SS later.
#' @param N a vector indicating the number of individuals to simulate
#' @param lam0 the detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the K x 2 matrix of trap locations
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @param tf an optional J x K matrix of dynamic trap operation indicating when and where 0, 1, or 2 cams were deployed
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from a camera trap SCR study. The extent of the state space is controlled by "buff", which buffers the
#' minimum and maximum X and Y extents.  Therefore, it probably only makes sense for square or rectangular grids.  Functionality
#' for user-provided polygonal state spaces will be added in the future.
#' @author Ben Augustine
#' @export
simOpenSPIM <-
  function(N=c(40,60,80),gamma=NULL,phi=rep(0.8),lam01=0.2,lam02=0.2,sigma=0.50,K=rep(10,3),X=X,t=3,M=M,sigma_t=NULL,buff=3,obstype="bernoulli",tf=NA){
    library(abind)
    #Check for user errors
    if(length(N)==1&is.null(gamma)){
      stop("Must provide gamma if length(N)==1")
    }
    if(length(N)==t&!is.null(gamma)){
      stop("Do not provide gamma if length(N)=t")
    }
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    if(length(X)!=t){
      stop("Must supply a X (trap locations) for each year")
    }
    storeparms=list(N=N,gamma=gamma,lam01=lam01,lam02=lam02,sigma=sigma,phi=phi)
    J=unlist(lapply(X,nrow))
    maxJ=max(J)
    maxK=max(K)
    #spim stuff
    for(i in 1:t){
      if(ncol(X[[i]])!=3){
        stop("X must have 3 columns for each year, X, Y, # cameras")
      }
    }
    if(!is.null(dim(tf))){
      if(!all(dim(tf)==c(J,K))){
        stop("tf must be of dimension J x K")
      }
      if(!all(tf==1|tf==2|tf==0)){
        stop("elements of tf must be 0, 1, or 2")
      }
      for(i in 1:t){
        twos=which(apply(tf[,,i],1,function(x){any(x==2)}))
        if(!all(twos%in%which(X[[i]][,3]==2))){
          stop("trap object X must have a 2 in 3rd column if tf indicates it had 2 cams operational at some point during survey")
        }
        if(any(twos%in%which(X[[i]][,3]==1))){
          stop("trap object must have a 1 in 3rd column if tf indicates it never had 2 cams operational")
        }
      }
      usetf=TRUE
    }else{
      usetf=FALSE
      tf=array(NA,dim=c(maxJ,maxK,t))
      for(l in 1:t){
        tf[1:J[l],,l]=X[[l]][,3]
      }
    }

    #Get state space extent such that buffer is at least buff on all grids
    #minmax=rbind(apply(X[[1]],2,min),apply(X[[1]],2,max))
    minmax=array(NA,dim=c(length(X),2,2))
    for(i in 1:length(X)){
      minmax[i,,]=rbind(apply(X[[i]][,1:2],2,min),apply(X[[i]][,1:2],2,max))
    }
    xylim=rbind(apply(minmax,3,min),apply(minmax,3,max))
    xlim=xylim[,1]
    ylim=xylim[,2]
    # # simulate a population of meta ACS and ACs, calculate D and lamd
    mu<- cbind(runif(M, xlim[1]-buff,xlim[2]+buff), runif(M,ylim[1]-buff,ylim[2]+buff))
    s=array(NA,dim=c(M,t,2))
    D=array(NA,dim=c(M,maxJ,t))
    lamd=array(NA,dim=c(M,maxJ,t,2))
    if(length(lam01)==1){
      lam01=rep(lam01,t)
    }
    if(length(lam02)==1){
      lam02=rep(lam02,t)
    }
    if(length(sigma)==1){
      sigma=rep(sigma,t)
    }
    if(is.null(sigma_t)){#no movement
      for(i in 1:t){
        s[,i,]=mu
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]][,1:2])
        lamd[,,i,1]=lam01[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        lamd[,,i,2]=lam02[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }
    }else{
      for(i in 1:t){#meta mu movement
        countout=0
        for(j in 1:M){
          out=1
          while(out==1){
            s[j,i,1]=rnorm(1,mu[j,1],sigma_t)
            s[j,i,2]=rnorm(1,mu[j,2],sigma_t)
            if(s[j,i,1] < xlim[2]+buff & s[j,i,1] > xlim[1]-buff & s[j,i,2] < ylim[2]+buff & s[j,i,2] > ylim[1]-buff){
              out=0
            }
            countout=countout+1
            if(countout==5000){#So you don't end up in an infinite loop if you set sigma_t too large relative to state space size
              stop("5000 ACs proposed outside of state space. Should probably reconsider settings.")
            }
          }
        }
        D[,1:nrow(X[[i]]),i]=e2dist(s[,i,],X[[i]])
        lamd[,,i,1]=lam01[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        lamd[,,i,2]=lam01[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }
    }
    psi=N[1]/M
    #Calculate (per capita) gamma or N
    if(length(phi)==1){
      phi=rep(phi,t-1)
    }

    if(is.null(gamma)){
      gamma=rep(NA,t-1)
      for(l in 2:t){
        gamma[l-1]=(N[l]-phi[t-1]*N[l-1])/(N[l-1])
      }
      storeparms$gamma=gamma
    }else{
      if(length(gamma)==1){
        gamma=rep(gamma,t-1)
      }
      N=c(N,rep(NA,t-1))
      for(l in 2:t){
        N[l]=round(N[l-1]*phi[t-1]+N[l-1]*gamma[l-1])
      }
    }
    #####Population Dynamics############# Stolen from Richard Chandler
    z=a=matrix(NA,M,t)
    # z[,1] <- rbinom(M, 1, EN1/M)
    z[1:N[1],1]=1
    z[(N[1]+1):M]=0
    a[,1]= 1-z[,1] # Available to be recruited?
    gamma.prime=rep(NA,4)
    for(i in 2:t) {
      ER <- sum(z[,i-1])*gamma[i-1] # Expected number of recruits
      A <- sum(a[,i-1]) # nAvailable to be recruited
      if(ER>A){
        stop("M is too low. There aren't any individuals left to be recruited")
      }
      gamma.prime[i-1] <- ER/A # individual-level recruitment *probability*
      Ez <- z[,i-1]*phi[i-1] + a[,i-1]*gamma.prime[i-1]
      z[,i] <- rbinom(M, 1, Ez)
      a[,i] <- apply(z[,1:i]<1, 1, all)
    }

    #######Capture process######################
    # Simulate encounter history
    left=right=both <-array(0,dim=c(M,maxJ,maxK,t))
    if(obstype=="bernoulli"){
      pd=cellprobsOpen(lamd)
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            for(k in 1:K[l]){
              left[i,j,k,l]=rbinom(1,1,(z[i,l]==1)*((tf[j,k,l]==1)*pd[i,j,l,1]+(tf[j,k,l]==2)*(2*pd[i,j,l,1]-pd[i,j,l,1]^2))) #single side p. two chances for capture with 2 cameras
              right[i,j,k,l]=rbinom(1,1,(z[i,l]==1)*((tf[j,k,l]==1)*pd[i,j,l,1]+(tf[j,k,l]==2)*(2*pd[i,j,l,1]-pd[i,j,l,1]^2))) #single side p
              both[i,j,k,l]=rbinom(1,1,(z[i,l]==1)*(pd[i,j,l,2])*(tf[j,k,l]==2)*1)  #both side lambda multiplied by indicator for 2 traps at site
            }
          }
        }
      }
    }else if(obstype=="poisson"){
      for(l in 1:t){
        for(i in 1:M){
          for(j in 1:J[l]){
            left[i,j,k,l]=rpois(1,(z[i,l]==1)*(tf[j,k,l]*lamd[i,j,l,1])) #single side lambda multiplied by number of traps at site
            right[i,j,k,l]=rpois(1,(z[i,l]==1)*(tf[j,k,l]*lamd[i,j,l,1])) #single side lambda multiplied by number of traps at site
            both[i,j,k,l]=rpois(1,(z[i,l]==1)*(lamd[i,j,l,2])*(tf[j,k,l]==2)*1)  #both side lambda multiplied by indicator for 2 traps at site          }
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    ######Process data#############
    IDknown=which(apply(both,1,sum)>0)
    Nknown=length(IDknown)
    n=rep(NA,t)
    nside=matrix(NA,ncol=3,nrow=t)
    nall=sum(apply(both+left+right,1,sum)>0)
    for(l in 1:t){
      n[l]=sum(apply(both[,,,l]+left[,,,l]+right[,,,1],1,sum)>0)
      nside[l,]=c(sum(apply(both[,,,l],1,sum)>0),sum(apply(left[,,,l],1,sum)>0),sum(apply(right[,,,l],1,sum)>0))
    }
    #Count true spatial recaps
    y2D=apply(both+left+right,c(1,2,4),sum)
    scaps=matrix(NA,ncol=t,nrow=M)
    nscap=sumscap=rep(NA,t)
    for(l in 1:t){
      scaps[,l]=rowSums(1*(y2D[,,l]>0))
      scaps[scaps[,l]>0,l]=scaps[scaps[,l]>0,l]-1
      nscap[l]=sum(scaps[,l]>0)
      sumscap[l]=sum(scaps[,l])
    }
    Srecap=data.frame(nscap=nscap,sumscap=sumscap)
    # relabel data sets if needed so that the left side always has more individuals than right after removing known IDs
    nl.tmp<- sum(rowSums(left)>0&rowSums(both)==0)
    nr.tmp<- sum(rowSums(right)>0&rowSums(both)==0)
    if(nr.tmp > nl.tmp){
      a<- left
      left<- right
      right<- a
    }
    # Move known individuals(after photo ID) to beginning
    if(Nknown>0&Nknown<M){
      IDunknown=which(is.na(match(1:M,IDknown)))
      newleft<- abind(left[IDknown,,,],left[IDunknown,,,],along=1)
      newright<- abind(right[IDknown,,,],right[IDunknown,,,],along=1)
      newboth<- abind(both[IDknown,,,],both[IDunknown,,,],along=1)
      s=s[c(IDknown,IDunknown),,]
      mu=mu[c(IDknown,IDunknown),]
      left<- newleft
      right<- newright
      both<- newboth
      #Recalculate after sorting
      IDknown=1:length(IDknown)
    }

    #Add all zero dimensions so both left and right can be added together
    # add=array(0,dim=c(M,maxJ,maxK,t))
    # both=abind(both,add,add,along=0)
    # left=abind(add,left,add,along=0)
    # right=abind(add,add,right,along=0)
    # both=aperm(both,c(2,1,3,4,5))
    # left=aperm(left,c(2,1,3,4,5))
    # right=aperm(right,c(2,1,3,4,5))

    #Remove all 0 capture histories.  Don't remove right and left all 0 histories for known IDs
    both<- both[IDknown,,,]
    if(length(dim(both))==3){#if there is 1 both guy need to keep 4d
      both=array(both,dim=c(1,dim(both)))
    }
    if(length(dim(left))==3){#if there is 1 left guy need to keep 4d
      left=array(left,dim=c(1,dim(left)))
    }
    if(length(dim(right))==3){#if there is 1 right guy need to keep 4d
      right=array(right,dim=c(1,dim(right)))
    }
    lcap=which(apply(left,1,sum)>0)
    rcap=which(apply(right,1,sum)>0)
    lorder=sort(unique(c(lcap,IDknown)))
    rorder=sort(unique(c(rcap,IDknown)))
    left<- left[lorder,,,]
    right<- right[rorder,,,]
    simID_L=lorder #keep these for sim plots
    simID_R=rorder

    #Renumber partial ID guys and put at end
    fix=unique(c(lorder[!lorder%in%IDknown],rorder[!rorder%in%IDknown]))
    if(length(fix)>0&Nknown>0){
      replace=length(IDknown)+1
      for(i in 1:length(fix)){
        if(fix[i]%in%lorder){
          lorder[which(lorder==fix[i])]=replace
        }
        if(fix[i]%in%rorder){
          rorder[which(rorder==fix[i])]=replace
        }
        replace=replace+1
      }

      ID_L=lorder
      ID_R=rorder
    }else{
      ID_L=ID_R=1:M
    }
    # #Build all possible observed data sets
    B3D=apply(both,c(1,2,4),sum)
    L3D=apply(left,c(1,2,4),sum)
    R3D=apply(right,c(1,2,4),sum)
    known=dim(both)[1]
    if(known>0){
      BLR=1*((both+left[1:known,,,]+right[1:known,,,])>0) #keep boths and lefts and rights for boths
      if(length(dim(BLR))==3){
        BLR=array(BLR,dim=c(1,dim(BLR)))
      }
      if(dim(left)[1]>known){
        BLRL=abind(BLR,left[(known+1):dim(left)[1],,,],along=1) #Add unknown lefts
      }else{
        BLRL=BLR
      }
      if(dim(right)[1]>known){
        BLRR=abind(BLR,right[(known+1):dim(right)[1],,,],along=1) #Add unknown rights
      }else{
        BLRR=BLR
      }
    }else{ #no known guys
      BLR=both
      if(dim(left)[1]>known){
        BLRL=left #Add unknown lefts
      }else{
        BLRL=BLR
      }
      if(dim(right)[1]>known){
        BLRR=right #Add unknown rights
      }else{
        BLRR=BLR
      }
      if(length(dim(BLRL))==3){
        BLRL=array(BLRL,dim=c(1,dim(BLRL)))
      }
      if(length(dim(BLRR))==3){
        BLRR=array(BLRR,dim=c(1,dim(BLRR)))
      }
    }
    BLR3D=apply(BLR,c(1,2,4),sum)
    BLRL3D=apply(BLRL,c(1,2,4),sum)
    BLRR3D=apply(BLRR,c(1,2,4),sum)
    y.obs=list(BLRL=BLRL,BLRR=BLRR)
    if(length(storeparms$gamma)==1){
      gamma=gamma[1]
    }else{
      gamma=gamma
    }
    if(is.null(sigma_t)){
      if(usetf==TRUE){
      out<-list(left=left,right=right,both=both,IDknown=IDknown,ID_L=ID_L,ID_R=ID_R,s=s,X=X,tf=tf,K=K,buff=buff,J=J,
                EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,Srecap=Srecap,simID_L=simID_L,simID_R=simID_R)
      }else{
        out<-list(left=left,right=right,both=both,IDknown=IDknown,ID_L=ID_L,ID_R=ID_R,s=s,X=X,K=K,buff=buff,J=J,
                  EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,Srecap=Srecap,simID_L=simID_L,simID_R=simID_R)
      }
    }else{
      if(usetf==TRUE){
        out<-list(left=left,right=right,both=both,IDknown=IDknown,ID_L=ID_L,ID_R=ID_R,mu=mu,s=s,X=X,tf=tf,K=K,buff=buff,J=J
                ,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,Srecap=Srecap,simID_L=simID_L,simID_R=simID_R)
      }else{
        out<-list(left=left,right=right,both=both,IDknown=IDknown,ID_L=ID_L,ID_R=ID_R,mu=mu,s=s,X=X,K=K,buff=buff,J=J
                  ,EN=N,N=colSums(z),z=z,gamma=gamma,phi=storeparms$phi,obstype=obstype,Srecap=Srecap,simID_L=simID_L,simID_R=simID_R)
      }
    }
    return(out)
  }
