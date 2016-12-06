SCRmcmcOpenSPIM <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200,swap=10,swap.tol=1, inits=inits,
           proppars=list(lam01=0.05,lam01=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,joint=FALSE){
    library(abind)
    t=dim(data$both)[5]
    both<-data$both
    left<-data$left
    right<-data$right
    if(dim(right)[1]>dim(left)[1]){
      storeL=left[,2,,,]
      storeR=right[,3,,,]
      dimL=dim(left)
      left=array(0,dim=dim(right))
      right=array(0,dim=dimL)
      left[,2,,,]=storeR
      right[,3,,,]=storeL
      warning("Right side data set larger than left so I switched them for convenience")
    }
    X<-data$X
    J<-data$J
    maxJ=max(J)
    K<-data$K
    maxK=max(K)
    IDknown=data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(left)[1]-Nfixed
    nright<-dim(right)[1]-Nfixed
    ID_L=data$ID_L
    ID_R=data$ID_R
    ####Open pop Error checks, add SPIM checks later
    if(length(K)!=t){
      stop("Must supply a K for each year")
    }
    metamu="sigma_t"%in%names(inits)
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
    lam01<- inits$lam01
    lam02<- inits$lam02
    sigma<- inits$sigma
    sigma_t=inits$sigma_t
    gamma=inits$gamma
    phi=inits$phi
    psi=inits$psi
    if(!length(lam01)%in%c(1,t)){
      stop("Input either 1 or t initial values for lam0")
    }
    if(!length(lam02)%in%c(1,t)){
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
    if(length(lam01)!=length(proppars$lam01)){
      stop("Must supply a tuning parameter for each lam0")
    }
    if(length(lam02)!=length(proppars$lam02)){
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

    #Figure out what needs to be updated
    uplam01=uplam02=upIDs=TRUE
    if(lam01==0){
      uplam01=FALSE
      upIDs=FALSE
    }
    all1=rep(FALSE,t)
    for(l in 1:t){
      if(all(X[[l]][,3]==1)){
        all1[l]=TRUE
      }
    }
    if(all(all1)|lam02==0){
      uplam02=FALSE
    }
    if(upIDs==TRUE&(nleft==0&nright==0)){  #not sure what will happen if we only have left only or right only
      upIDs=FALSE
    }

    #sort to minimize distance between initial matches. Skip if no single sides.
    if(nleft>0|nright>0){
      IDs<- LRmatchOpen(M=M,left=left, nleft=nleft, right=right, nright=nright, X, Nfixed=Nfixed)
      #Add unused augmented indivuals back in
      notusedL<- (1:M)[is.na(match(1:M,IDs$ID_L))]
      ID_L<-c(IDs$ID_L,notusedL)
      notusedR<- (1:M)[is.na(match(1:M,IDs$ID_R))]
      ID_R<-c(IDs$ID_R,notusedR)
    }else{
      ID_R=ID_L=1:M
    }
    #augment data
    both<- abind(both,array(0, dim=c( M-dim(both)[1],3,maxJ,maxK, t)), along=1)
    left<- abind(left,array(0, dim=c( M-dim(left)[1],3,maxJ,maxK, t)), along=1)
    right<- abind(right,array(0, dim=c( M-dim(right)[1],3,maxJ,maxK, t)), along=1)

    #Make initial complete data set
    tmpdata<- both + left[order(ID_L),,,,] + right[order(ID_R),,,,]
    tmpdata<- apply(tmpdata,c(1,3,5),sum)
    # z=1*(apply(tmpdata,1,sum)>0)
    # z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )

    #Initialize z, r, and a consistent with y
    known.matrix=1*(apply(tmpdata,c(1,3),sum)>0)#make z consistent with y
    if(t>2){
      for(l in 2:(t-1)){#Turn on zeros with 1's on either side
        known.matrix[known.matrix[,l]==0&known.matrix[,l-1]==1&rowSums(matrix(known.matrix[,(l+1):t],nrow=M))>0,l]=1
      }
    }
    z=known.matrix
    # r=array(0,dim=dim(z))
    #turn on z's with caps on either side so we know they were in pop
    knownguys=which(rowSums(tmpdata)>0)
    unknownguys=setdiff(1:M,knownguys)
    for(i in knownguys){
      if(sum(z[i,]>0)){
        idx=which(z[i,]==1)
        z[i,min(idx):max(idx)]=1
      }
    }
    #turn on some augmented guys for z1
    # z[(n+1):M,1]=rbinom(M-n,1,psi)#add augmented guys to t=1 with psi.. not enough
    z1deal=psi*M-sum(z[,1])
    if(z1deal<0){
      stop("initial psi is too small given M to turn on any uncaptured z[,1] guys ")
    }
    z[sample(unknownguys,z1deal),1]=1
    a=matrix(1,nrow=M,ncol=t) #a is available to be recruited
    a[which(z[,1]==1),]=0#turn off guys caught year 1
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
      nrecruit=round(sum(z[,i-1])*gammasim[i-1])#how many should we recruit?
      recruits=which(z[,i-1]==0&z[,i]==1)#who is already recruited based on y constraints?
      a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleft=nrecruit-length(recruits)
      if(nleft>0){
        cands=setdiff(which(a[,i-1]==1),recruits)#Who is available to recruit?
        cands=cands[!cands%in%knownguys]#Don't mess with known guys
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        pick=sample(cands,nleft)
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Can't initialize all E[recruits] given initial value for gamma. Should probably raise M?")
      }
      #Survival
      nlive=rbinom(1,sum(z[,i-1]),phisim[i-1])#How many to live
      livers=which(z[,i-1]==1&z[,i]==1)
      nleft=nlive-length(livers)
      a[livers,i]=0
      if(nleft>0){
        cands=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1)
        cands=cands[!cands%in%knownguys]
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
    gamma.prime.cand=gamma.prime
    ll.z.cand=ll.z
    s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(tmpdata)>0) #switch for those actually caught
    for(i in idx){
      trps=matrix(0,nrow=0,ncol=2)
      for(l in 1:t){ #loop over t to get all cap locs
        trps<- rbind(trps,X[[l]][tmpdata[i,,l]>0,1:2])
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
        idx=which(rowSums(tmpdata[,,l])>0) #switch for those actually caught
        for(i in 1:M){
          if(i%in%idx){
            trps<- X[[l]][tmpdata[i,,l]>0,1:2]
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
      #Count z==0
      ll.s2=(dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE))
    }
    # some objects to hold the MCMC simulation output
    if(niter<(nburn)){
      stop("niter is smaller than nburn")
    }
    nstore=(niter-nburn)/nthin
    if((nburn)%%nthin!=0){
      nstore=nstore+1
    }
    if(length(lam01)==t){
      lam01names=paste("lam01",1:t,sep="")
    }else{
      lam01names="lam01"
    }
    if(length(lam02)==t){
      lam02names=paste("lam02",1:t,sep="")
    }else{
      lam02names="lam02"
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
      out<-matrix(NA,nrow=nstore,ncol=length(lam01)+length(lam02)+length(sigma)+length(gamma)+length(phi)+t+1)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"sigma_t")
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else{
      out<-matrix(NA,nrow=nstore,ncol=length(lam01)+length(lam02)+length(sigma)+length(gamma)+length(phi)+t)
      colnames(out)<-c(lam01names,lam02names,sigmanames,gammanames,phinames,Nnames)
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
    }
    ID_Lout<- ID_Rout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    D=lamd1=lamd2=lik.curr=array(NA,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(i in 1:t){
      D[,1:nrow(X[[i]]),i]=e2dist(s2[,i,],X[[i]])
      if(length(lam01)==t&length(sigma)==t){
        lamd1[,,i]=lam01[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        lamd2[,,i]=lam02[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam01)==1&length(sigma)==t){
        lamd1[,,i]=lam01*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
        lamd2[,,i]=lam02*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam01)==t&length(sigma)==1){
        lamd1[,,i]=lam01[i]*exp(-D[,,i]^2/(2*sigma*sigma))
        lamd2[,,i]=lam02[i]*exp(-D[,,i]^2/(2*sigma*sigma))
      }else{
        lamd1[,,i]=lam01*exp(-D[,,i]^2/(2*sigma*sigma))
        lamd2[,,i]=lam02*exp(-D[,,i]^2/(2*sigma*sigma))
      }
    }

    y.both<- apply(both[,1,,,], c(1,2,4), sum)
    y.left<- apply(left[order(ID_L),2,,,], c(1,2,4), sum)
    y.right<- apply(right[order(ID_R),3,,,], c(1,2,4), sum)
    zero.guys<- apply(y.both+y.left + y.right ,c(1,3),sum) == 0

    #Calculate ll for observation model
    pd1=1-exp(-lamd1)
    pd2=1-exp(-lamd2)
    trapno=ones=twos=array(NA,dim=c(M,maxJ,t))
    for(l in 1:t){
      trapno[,,l]=matrix(rep(X[[l]][,3],M),nrow=M,byrow=TRUE) #trap number multiplier for left and right captures
    }
    ones=trapno==1
    twos=trapno==2
    ll.y.both=ll.y.left=ll.y.right=array(NA,dim=dim(y.both))
    for(l in 1:t){
      ll.y.both[,,l]=dbinom(y.both[,,l],K[l],pd2[,,l]*z[,l],log=TRUE)
      ll.y.both[,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
      ll.y.left[,,l]=dbinom(y.left[,,l],K[l],z[,l]*(ones[,,l]*pd1[,,l]+twos[,,l]*(2*pd1[,,l]-pd1[,,l]*pd1[,,l])),log=TRUE)
      ll.y.right[,,l]=dbinom(y.right[,,l],K[l],z[,l]*(ones[,,l]*pd1[,,l]+twos[,,l]*(2*pd1[,,l]-pd1[,,l]*pd1[,,l])),log=TRUE)
    }
    ll.y.both.cand=ll.y.both
    ll.y.left.cand=ll.y.left
    ll.y.right.cand=ll.y.right
    ll.y.both.t.sum=ll.y.both.cand.t.sum=apply(ll.y.both,3,sum) #ll summed for each year
    ll.y.left.t.sum=ll.y.left.cand.t.sum=apply(ll.y.left,3,sum) #ll summed for each year
    ll.y.right.t.sum=ll.y.right.cand.t.sum=apply(ll.y.right,3,sum) #ll summed for each year
    ll.y.both.sum=sum(ll.y.both.t.sum) #full ll sum
    ll.y.left.sum=sum(ll.y.left.t.sum) #full ll sum
    ll.y.right.sum=sum(ll.y.right.t.sum) #full ll sum
    lamd1.cand=lamd1
    lamd2.cand=lamd2
    pd1.cand=pd1
    pd2.cand=pd2

    #Check ll.y.sum for inf
    if(!is.finite(ll.y.both.sum)){
      stop("Both side detection function starting values produce -Inf log likelihood values. Try increasing sigma and/or lam02")
    }
    if(!is.finite(ll.y.left.sum)|!is.finite(ll.y.left.sum)){
      stop("Single side detection function starting values produce -Inf log likelihood values. Try increasing sigma and/or lam01")
    }
    if(joint==TRUE){
      #Figure out all possible z histories
      zpossible=cbind(c(1,1,0),c(1,0,1))
      if (t > 2) {
        for (i in (3:t)) {
          zpossible=cbind(c(rep(1,2^(i-1)),rep(0,((2^(i-1))-1))), rbind(zpossible, rep(0, (i-1)), zpossible))
        }
      }
      #remove zombie histories
      illegal=rep(FALSE,nrow(zpossible))
      for(l in 1:(t-2)){
        latecaps=rep(0,nrow(zpossible))
        for(l2 in (l+2):(t)){
          latecaps=latecaps+zpossible[,l2]
        }
        illegal=illegal|(zpossible[,l]==1&zpossible[,l+1]==0&latecaps>0)
      }
      zpossible=zpossible[which(illegal==FALSE),]
      nzpossible=nrow(zpossible)
      apossible=matrix(1,nrow=nzpossible,ncol=t)
      apossible[zpossible[,1]==1,1]=0
      for(l in 2:t){
        apossible[apossible[,l-1]==1&zpossible[,l]==1,l]=0
        apossible[apossible[,l-1]==0,l]=0
      }
      Ezpossible=matrix(NA,nrow=nzpossible,ncol=t-1)
      ll_z_possible=matrix(NA,nrow=nzpossible,ncol=t)
    }

    for(iter in 1:niter){
      ll.y.both.t.sum=apply(ll.y.both,3,sum) #only needed for detection parameters, changes in z and AC updates
      ll.y.left.t.sum=apply(ll.y.left,3,sum)
      ll.y.right.t.sum=apply(ll.y.right,3,sum)
      ll.y.both.sum=sum(ll.y.both.t.sum)
      ll.y.left.sum=sum(ll.y.left.t.sum)
      ll.y.right.sum=sum(ll.y.right.t.sum)
      # Update lam01
      if(uplam01){
        if(length(lam01)==t){ #if lam0 is year-specific
          for(l in 1:t){
            lam01.cand<- rnorm(1,lam01[l],proppars$lam01[l])
            if(lam01.cand > 0){
              if(length(sigma)==t){#if sigma is year specific
                lamd1.cand[,,l]<- lam01.cand*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
              }else{#fixed sigma
                lamd1.cand[,,l]<- lam01.cand*exp(-D[,,l]^2/(2*sigma*sigma))
              }
              pd1.cand[,,l]=1-exp(-lamd1.cand[,,l])
              ll.y.left.cand[,,l]=dbinom(y.left[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
              ll.y.right.cand[,,l]=dbinom(y.right[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
              ll.y.left.cand.t.sum[l]=sum(ll.y.left.cand[,,l])
              ll.y.right.cand.t.sum[l]=sum(ll.y.right.cand[,,l])
              if(runif(1) < exp((ll.y.left.cand.t.sum[l]+ll.y.right.cand.t.sum[l]) -(ll.y.right.t.sum[l]+ll.y.right.t.sum[l]))){
                lam01[l]<- lam01.cand
                lamd1[,,l]=lamd1.cand[,,l]
                pd1[,,l]=pd1.cand[,,l]
                ll.y.left[,,l]=ll.y.left.cand[,,l]
                ll.y.right[,,l]=ll.y.right.cand[,,l]
                ll.y.left.t.sum[l]=ll.y.left.cand.t.sum[l]
                ll.y.right.t.sum[l]=ll.y.right.cand.t.sum[l]
              }
            }
          }
          ll.y.left.sum=sum(ll.y.left.t.sum)
          ll.y.right.sum=sum(ll.y.right.t.sum)
        }else{#fixed lam01
          lam01.cand<- rnorm(1,lam01,proppars$lam01)
          if(lam01.cand > 0){
            if(length(sigma)==t){#if sigma is year specific
              for(l in 1:t){
                lamd1.cand[,,l]<- lam01.cand*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
              }
            }else{#fixed sigma
              lamd1.cand<- lam01.cand*exp(-D^2/(2*sigma*sigma))
            }
            pd1.cand=1-exp(-lamd1.cand)
            for(l in 1:t){
              ll.y.left.cand[,,l]=dbinom(y.left[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
              ll.y.right.cand[,,l]=dbinom(y.right[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
            }
            ll.y.left.cand.t.sum=apply(ll.y.left.cand,3,sum)
            ll.y.right.cand.t.sum=apply(ll.y.right.cand,3,sum)
            ll.y.left.cand.sum=sum(ll.y.left.cand.t.sum)
            ll.y.right.cand.sum=sum(ll.y.right.cand.t.sum)
            if(runif(1) < exp((ll.y.left.cand.sum+ll.y.right.cand.sum )- (ll.y.left.sum+ll.y.right.sum))){
              lam01= lam01.cand
              lamd1=lamd1.cand
              pd1=pd1.cand
              ll.y.left=ll.y.left.cand
              ll.y.right=ll.y.right.cand
              ll.y.left.sum=ll.y.left.cand.sum
              ll.y.right.sum=ll.y.right.cand.sum
              ll.y.left.t.sum=ll.y.left.cand.t.sum
              ll.y.right.t.sum=ll.y.right.cand.t.sum
            }
          }
        }
      }
      # Update lam02
      if(uplam02){
        if(length(lam02)==t){ #if lam0 is year-specific
          for(l in 1:t){
            lam02.cand<- rnorm(1,lam02[l],proppars$lam02[l])
            if(lam02.cand > 0){
              if(length(sigma)==t){#if sigma is year specific
                lamd2.cand[,,l]<- lam02.cand*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
              }else{#fixed sigma
                lamd2.cand[,,l]<- lam02.cand*exp(-D[,,l]^2/(2*sigma*sigma))
              }
              pd2.cand[,,l]=1-exp(-lamd2.cand[,,l])
              ll.y.both.cand[,,l]=dbinom(y.both[,,l],K[l],pd2.cand[,,l]*z[,l],log=TRUE)
              ll.y.both.cand[,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
              ll.y.both.cand.t.sum[l]=sum(ll.y.both.cand[,,l])
              if(runif(1) < exp((ll.y.both.cand.t.sum[l]) -(ll.y.both.t.sum[l]))){
                lam02[l]<- lam02.cand
                lamd2[,,l]=lamd2.cand[,,l]
                pd2[,,l]=pd2.cand[,,l]
                ll.y.both[,,l]=ll.y.both.cand[,,l]
                ll.y.both.t.sum[l]=ll.y.both.cand.t.sum[l]
              }
            }
          }
          ll.y.both.sum=sum(ll.y.both.t.sum)

        }else{#fixed lam02
          lam02.cand<- rnorm(1,lam02,proppars$lam02)
          if(lam02.cand > 0){
            if(length(sigma)==t){#if sigma is year specific
              for(l in 1:t){
                lamd2.cand[,,l]<- lam02.cand*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
              }
            }else{#fixed sigma
              lamd2.cand<- lam02.cand*exp(-D^2/(2*sigma*sigma))
            }
            pd2.cand=1-exp(-lamd2.cand)
            for(l in 1:t){
              ll.y.both.cand[,,l]=dbinom(y.both[,,l],K[l],pd2.cand[,,l]*z[,l],log=TRUE)
              ll.y.both.cand[,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
            }
            ll.y.both.cand.t.sum=apply(ll.y.both.cand,3,sum)
            ll.y.both.cand.sum=sum(ll.y.both.cand.t.sum)
            if(runif(1) < exp((ll.y.both.cand.sum )- (ll.y.both.sum))){
              lam02= lam02.cand
              lamd2=lamd2.cand
              pd2=pd2.cand
              ll.y.both.t.sum=ll.y.both.cand.t.sum
              ll.y.both.sum=ll.y.both.cand.sum
            }
          }
        }
      }
      #Update sigma
      if(length(sigma)==t){ #if sigma is year-specific
        for(l in 1:t){
          sigma.cand<- rnorm(1,sigma[l],proppars$sigma[l])
          if(sigma.cand > 0){
            if(length(lam01)==t){#if lam01 and lam02 are year specific
              lamd1.cand[,,l]<- lam01[l]*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
              lamd2.cand[,,l]<- lam02[l]*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
            }else{#fixed lam0
              lamd1.cand[,,l]<- lam01*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
              lamd2.cand[,,l]<- lam02*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
            }
            #single sides
            pd1.cand[,,l]=1-exp(-lamd1.cand[,,l])
            ll.y.left.cand[,,l]=dbinom(y.left[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
            ll.y.right.cand[,,l]=dbinom(y.right[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
            ll.y.left.cand.t.sum[l]=sum(ll.y.left.cand[,,l])
            ll.y.right.cand.t.sum[l]=sum(ll.y.right.cand[,,l])

            #both sides
            pd2.cand[,,l]=1-exp(-lamd2.cand[,,l])
            ll.y.both.cand[,,l]=dbinom(y.both[,,l],K[l],pd2.cand[,,l]*z[,l],log=TRUE)
            ll.y.both.cand[,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
            ll.y.both.cand.t.sum[l]=sum(ll.y.both.cand[,,l])
            ll.y.cand.t.sum[l]=sum(ll.y.cand[,,l])#just 1 year

            if(runif(1) < exp((ll.y.both.cand.t.sum[l]+ll.y.left.cand.t.sum[l]+ll.y.right.cand.t.sum[l]) -(ll.y.both.t.sum[l]+ll.y.left.t.sum[l]+ll.y.right.t.sum[l]))){
              sigma[l]<- sigma.cand
              lamd1[,,l]=lamd1.cand[,,l]
              pd1[,,l]=pd1.cand[,,l]
              ll.y.left[,,l]=ll.y.left.cand[,,l]
              ll.y.right[,,l]=ll.y.right.cand[,,l]
              ll.y.left.t.sum[l]=ll.y.left.cand.t.sum[l]
              ll.y.right.t.sum[l]=ll.y.right.cand.t.sum[l]
              lamd2[,,l]=lamd2.cand[,,l]
              pd2[,,l]=pd2.cand[,,l]
              ll.y.both[,,l]=ll.y.both.cand[,,l]
              ll.y.both.t.sum[l]=ll.y.both.cand.t.sum[l]
            }
          }
        }
        ll.y.sum=sum(ll.y.t.sum)
      }else{#fixed sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          if(length(lam01)==t){#if lam01 and lam02 are year specific
            for(l in 1:t){
              lamd1.cand[,,l]<- lam01[l]*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
              lamd2.cand[,,l]<- lam02[l]*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
            }
          }else{#fixed lam0
            lamd1.cand<- lam01*exp(-D^2/(2*sigma.cand*sigma.cand))
            lamd2.cand<- lam02*exp(-D^2/(2*sigma.cand*sigma.cand))
          }
          pd1.cand=1-exp(-lamd1.cand)
          pd2.cand=1-exp(-lamd2.cand)
          for(l in 1:t){
            ll.y.left.cand[,,l]=dbinom(y.left[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
            ll.y.right.cand[,,l]=dbinom(y.right[,,l],K[l],z[,l]*(ones[,,l]*pd1.cand[,,l]+twos[,,l]*(2*pd1.cand[,,l]-pd1.cand[,,l]*pd1.cand[,,l])),log=TRUE)
            ll.y.both.cand[,,l]=dbinom(y.both[,,l],K[l],pd2.cand[,,l]*z[,l],log=TRUE)
            ll.y.both.cand[,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
          }
          ll.y.left.cand.t.sum=apply(ll.y.left.cand,3,sum)
          ll.y.right.cand.t.sum=apply(ll.y.right.cand,3,sum)
          ll.y.left.cand.sum=sum(ll.y.left.cand.t.sum)
          ll.y.right.cand.sum=sum(ll.y.right.cand.t.sum)
          ll.y.both.cand.t.sum=apply(ll.y.both.cand,3,sum)
          ll.y.both.cand.sum=sum(ll.y.both.cand.t.sum)
          if(runif(1) < exp((ll.y.both.cand.sum+ll.y.left.cand.sum+ll.y.right.cand.sum) - (ll.y.both.sum+ll.y.left.sum+ll.y.right.sum))){
            sigma<- sigma.cand
            lamd1=lamd1.cand
            lamd2=lamd2.cand
            pd1=pd1.cand
            pd2=pd2.cand
            ll.y.left=ll.y.left.cand
            ll.y.right=ll.y.right.cand
            ll.y.left.sum=ll.y.left.cand.sum
            ll.y.right.sum=ll.y.right.cand.sum
            ll.y.left.t.sum=ll.y.left.cand.t.sum
            ll.y.right.t.sum=ll.y.right.cand.t.sum
            ll.y.both.t.sum=ll.y.both.cand.t.sum
            ll.y.both.sum=ll.y.both.cand.sum
          }
        }
      }

      ## Update latent ID variables
      ## Candidate: swap IDs of one guy with some other guy. The two guys to be swapped are
      ## chosen randomly from the z=1 guys
      if(any(!zero.guys&z==0)){
        cat("some z = 0! going into left update",fill=TRUE)
        browser()
      }
      zAny=rowSums(z)>0 #zAny status cannot change while updating, just which z's are on
      #Swap left guys first
      if(nleft>0){
        # User inputs how many swaps to make on each iteration my specifying "swap"
        #map lefts to boths
        #candmap used to remove disallowed candidates.NA any lefts that map back to z=0 indices and boths that are z=0 indices
        map=cbind(1:M,ID_L)
        candmap=map
        candmap[1:Nfixed,]=NA #Don't swap out IDknown guys
        candmap[which(!zAny),1]=NA #Don't swap in z=0 guys.
        # candmap[which(!zAny),]=NA #Don't swap in z=0 guys.
        candmap[candmap[,2]%in%which(!zAny),2]=NA #Don't choose a guy1 that will make you swap in a z=0 guy
        #These are the guys that can come out and go in
        OUTcands=which(!is.na(candmap[,2]))
        INcands=which(!is.na(candmap[,1]))
        for(up in 1:swap){
          ### In this code here "guy1" and "guy2" are indexing "left-side encounter histories" whereas
          ### "s.swap.in" and "s.swap.out" are indexing (both-side guys, s, z)
          guy1<- sample(OUTcands, 1)
          s.swap.out<- map[guy1,2]
          # to find candidates for swapping look in vicinity of it's current both side membership
          dv<-  sqrt( (s1[s.swap.out,1]- s1[INcands,1])^2 + (s1[s.swap.out,2] - s1[INcands,2])^2 )
          ## This is a list of both side activity centers that could be assinged to the right side guy
          ## under consideration based on how far away they are.  swap.tol is the spatial tolerance
          # if no one around, code swaps guy 1 with guy 1 (no swap, but does calculations)
          possible<- INcands[dv < swap.tol]

          if(joint==FALSE){
            if(length(possible)==0){
              next
            }else{
              keep=rep(NA,length(possible))
              for(j in 1:length(possible)){
                keep[j]=all(z[s.swap.out,]==z[possible[j],])
              }
            }
            possible=possible[keep]
            # this is a particular value of s to swap for guy1.id
            if(length(possible)>1){
              s.swap.in <-  sample( possible, 1)
            }
            if(length(possible) ==1){
              s.swap.in <- possible
            }
            if(length(possible)==0) next #no guys close enough to swap
            jump.probability<- 1/length(possible)  # h(theta*|theta)
            # z[c(s.swap.in,s.swap.out),]


            #compute  h(theta|theta*)
            trash<-   sqrt( (s1[s.swap.in,1]- s1[INcands,1])^2 + (s1[s.swap.in,2] - s1[INcands,2])^2  )
            trash<-  INcands[trash < swap.tol]
            jump.back<-  1/length(trash)

            ##  Which left encounter history is currently associated with both guy s.swap.in?
            guy2<- which(map[,2]==s.swap.in)
            if(guy1==guy2) next
            newID<-ID_L
            newID[guy1]<- s.swap.in
            newID[guy2]<- s.swap.out

            ## recompute 'data' and compute likelihood components for swapped guys only
            y.left.tmp<- apply( left[order(newID),2,,,], c(1,2,4), sum)#Only y.left changed
            swapped=c(s.swap.out,s.swap.in)
            for(l in 1:t){
              ll.y.left.cand[swapped,,l]=dbinom(y.left.tmp[swapped,,l],K[l],z[swapped,l]*(ones[swapped,,l]*pd1[swapped,,l]+twos[swapped,,l]*(2*pd1[swapped,,l]-pd1[swapped,,l]*pd1[swapped,,l])),log=TRUE)
            }
            if(!is.finite(sum(ll.y.left.cand))){
              stop("ll.y.left not finite. Maybe swap.tol too large")
            }
            llswap.curr<- sum(ll.y.left[swapped,,l])
            llswap.cand<- sum(ll.y.left.cand[swapped,,l])
            lldiff=llswap.cand-llswap.curr
            #MH step
            if(runif(1)<exp(lldiff)*(jump.back/jump.probability) ){
              # lik.curr=lik.curr+lldiff #update likelihood
              y.left[swapped,,] <- y.left.tmp[swapped,,] #update left data
              ll.y.left[swapped,,]=ll.y.left.cand[swapped,,]
              ID_L<-newID #update data order
              map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
              zero.guys<- apply(y.both+y.left + y.right ,c(1,3),sum) == 0
              if(any(((colSums(y.both[guy1,,] + y.left.tmp[guy1,,]+ y.right[guy1,,])>0)==TRUE)&((z[guy1,]==0)==TRUE))){
                cat("error on guy1",fill=TRUE)
                browser()
              }
              if( any(((colSums(y.both[guy2,,] + y.left.tmp[guy2,,]+ y.right[guy2,,])>0)==TRUE)&((z[guy2,]==0)==TRUE))){
                cat("error on guy2",fill=TRUE)
                browser()
              }
            }
          }else{
            #joint update
            if(length(possible)>1){
              s.swap.in <-  sample( possible, 1)
            }
            if(length(possible) ==1){
              s.swap.in <- possible
            }
            if(length(possible)==0) next #no guys close enough to swap
            jump.probability<- 1/length(possible)  # h(theta*|theta)

            #compute  h(theta|theta*)
            trash<-   sqrt( (s1[s.swap.in,1]- s1[INcands,1])^2 + (s1[s.swap.in,2] - s1[INcands,2])^2  )
            trash<-  INcands[trash < swap.tol]
            jump.back<-  1/length(trash)

            ##  Which left encounter history is currently associated with both guy s.swap.in?
            guy2<- which(map[,2]==s.swap.in)
            if(guy1==guy2)next
            newID<-ID_L
            newID[guy1]<- s.swap.in
            newID[guy2]<- s.swap.out

            ## recompute 'data' and compute likelihood components for swapped guys only
            y.left.tmp<- apply( left[order(newID),2,,,], c(1,2,4), sum)#Only y.left changed
            swapped=c(s.swap.in,s.swap.out)

            ##update z before ll.y.left
            #Get likelihood for all possible z histories
            ll_z_possible[,1]=dbinom(zpossible[,1], 1, psi,log=TRUE)
            for(l in 2:t){
              Ezpossible[,l-1]=zpossible[,l-1]*phiuse[l-1] + apossible[,l-1]*gamma.prime[l-1]
              ll_z_possible[,l]=dbinom(zpossible[,l], 1, Ezpossible[,l-1],log=TRUE)
            }

            #Proposing y.left.tmp changes known.matrix
            swappedR=c(which(ID_R==swapped[1]),which(ID_R==swapped[2]))
            tmpdata.prop=both[swapped,,,,] + left[c(guy1,guy2),,,,] + right[swappedR,,,,]
            tmpdata.prop<- apply(tmpdata.prop,c(1,3,5),sum)
            known.matrix.prop=1*(apply(tmpdata.prop,c(1,3),sum)>0)

            #pick a new z based on updated known.matrix constraints
            zchoose=zcandguy=prop.prob=back.prob=rep(NA,2)
            zprop=matrix(NA,nrow=2,ncol=t)
            for(i2 in 1:2){
              #Zero out known matrix years for both swapped
              fixed=which(known.matrix.prop[i2,]==1)
              cancel=rep(1,nzpossible)
              for(i3 in 1:nzpossible){
                if(!all(fixed%in%which(zpossible[i3,]==1))){
                  cancel[i3]=0
                }
              }
              #new z stuff
              propto1=rowSums(exp(ll_z_possible))*(1*(cancel==1))
              propto=propto1/sum(propto1)
              zchoose[i2]=sample(1:nzpossible,1,prob=propto)
              zprop[i2,]=zpossible[zchoose[i2],]
              prop.prob[i2]=propto[zchoose[i2]]
              #old z stuff
              zcandguy[i2]=which(apply(zpossible,1,function(x){all(x==z[swapped[i2],])}))
              back.prob[i2]=propto[zcandguy[i2]]
            }

            #Because a and z changes, must update gamma.prime and Ez
            #Don't need to update all years every time, but not figuring that out for now
            aprop=apossible[zchoose,]

            #Calculate gamma.prime, Ez, and ll.z candidates
            ll.z.cand[swapped,1] <- dbinom(zprop[,1], 1, psi, log=TRUE)
            z.cand=z
            z.cand[swapped,]=zprop
            a.cand=a
            a.cand[swapped,]=aprop
            Ntmp=colSums(z.cand)
            for(l in 2:t){
              gamma.prime.cand[l-1]=(Ntmp[l-1]*gammause[l-1]) / sum(a.cand[,l-1])
              if(gamma.prime.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
                warning("Rejected z due to low M")
                next
              }
              Ez.cand[,l-1]=z.cand[,l-1]*phiuse[l-1] + a.cand[,l-1]*gamma.prime.cand[l-1]
              ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
            }

            #update new ll.y.left for new y.left and z
            for(l in 1:t){
              ll.y.left.cand[swapped,,l]=dbinom(y.left.tmp[swapped,,l],K[l],zprop[,l]*(ones[swapped,,l]*pd1[swapped,,l]+twos[swapped,,l]*(2*pd1[swapped,,l]-pd1[swapped,,l]*pd1[swapped,,l])),log=TRUE)
            }
            if(!is.finite(sum(ll.y.left.cand))){
              stop("ll.y.left not finite. Maybe swap.tol too large")
            }
            llyswap.curr<- sum(ll.y.left[swapped,,])
            llyswap.cand<- sum(ll.y.left.cand[swapped,,])
            llzswap.curr=sum(ll.z[swapped,1])+sum(ll.z[,2:t])
            llzswap.cand=sum(ll.z.cand[swapped,1])+sum(ll.z.cand[,2:t])
            #MH step
            if(runif(1)<exp((llyswap.cand+llzswap.cand)-(llyswap.curr+llzswap.curr))*((jump.back*prod(back.prob))/(jump.probability*prod(prop.prob)))){
              y.left[swapped,,] <- y.left.tmp[swapped,,] #update left data
              ll.y.left[swapped,,]=ll.y.left.cand[swapped,,]
              ID_L<-newID #update data order
              map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
              known.matrix[swapped,]=known.matrix.prop
              gamma.prime=gamma.prime.cand
              N=Ntmp
              Ez=Ez.cand
              ll.z[swapped,1]=ll.z.cand[swapped,1]
              ll.z[,2:t]=ll.z.cand[,2:t]
              zero.guys<- apply(y.both+y.left + y.right ,c(1,3),sum) == 0
              if(any(((colSums(y.both[guy1,,] + y.left.tmp[guy1,,]+ y.right[guy1,,])>0)==TRUE)&((z[guy1,]==0)==TRUE))){
                cat("error on guy1",fill=TRUE)
                browser()
              }
              if( any(((colSums(y.both[guy2,,] + y.left.tmp[guy2,,]+ y.right[guy2,,])>0)==TRUE)&((z[guy2,]==0)==TRUE))){
                cat("error on guy2",fill=TRUE)
                browser()
              }
            }
          }
        }
        #Do any captured guys have z==0?
        if(any(!zero.guys&z==0)){
          cat("coming out of ID_L update: some z = 0!",fill=TRUE)
          browser()
        }
      }

      #Repeat for rights
      if(nright>0){
        # User inputs how many swaps to make on each iteration my specifying "swap"
        #map lefts to boths
        #candmap used to remove disallowed candidates.NA any lefts that map back to z=0 indices and boths that are z=0 indices
        map=cbind(1:M,ID_R)
        candmap=map
        candmap[1:Nfixed,]=NA #Don't swap out IDknown guys
        candmap[which(!zAny),1]=NA #Don't swap in z=0 guys.
        # candmap[which(!zAny),]=NA #Don't swap in z=0 guys.
        candmap[candmap[,2]%in%which(!zAny),2]=NA #Don't choose a guy1 that will make you swap in a z=0 guy
        #These are the guys that can come out and go in
        OUTcands=which(!is.na(candmap[,2]))
        INcands=which(!is.na(candmap[,1]))
        for(up in 1:swap){
          ### In this code here "guy1" and "guy2" are indexing "left-side encounter histories" whereas
          ### "s.swap.in" and "s.swap.out" are indexing (both-side guys, s, z)
          guy1<- sample(OUTcands, 1)
          s.swap.out<- map[guy1,2]
          # to find candidates for swapping look in vicinity of it's current both side membership
          dv<-  sqrt( (s1[s.swap.out,1]- s1[INcands,1])^2 + (s1[s.swap.out,2] - s1[INcands,2])^2 )
          ## This is a list of both side activity centers that could be assinged to the right side guy
          ## under consideration based on how far away they are.  swap.tol is the spatial tolerance
          # if no one around, code swaps guy 1 with guy 1 (no swap, but does calculations)
          possible<- INcands[dv < swap.tol]
          if(joint==FALSE){
            if(length(possible)==0){
              next
            }else{
              keep=rep(NA,length(possible))
              for(j in 1:length(possible)){
                keep[j]=all(z[s.swap.out,]==z[possible[j],])
              }
            }
            possible=possible[keep]
            # this is a particular value of s to swap for guy1.id
            if(length(possible)>1){
              s.swap.in <-  sample( possible, 1)
            }
            if(length(possible) ==1){
              s.swap.in <- possible
            }
            if(length(possible)==0) next #no guys close enough to swap
            jump.probability<- 1/length(possible)  # h(theta*|theta)

            #compute  h(theta|theta*)
            trash<-   sqrt( (s1[s.swap.in,1]- s1[INcands,1])^2 + (s1[s.swap.in,2] - s1[INcands,2])^2  )
            trash<-  INcands[trash < swap.tol]
            jump.back<-  1/length(trash)

            ##  Which left encounter history is currently associated with both guy s.swap.in?
            guy2<- which(map[,2]==s.swap.in)
            if(guy1==guy2) next
            newID<-ID_R
            newID[guy1]<- s.swap.in
            newID[guy2]<- s.swap.out

            ## recompute 'data' and compute likelihood components for swapped guys only
            y.right.tmp<- apply( left[order(newID),3,,,], c(1,2,4), sum)#Only y.left changed
            swapped=c(s.swap.out,s.swap.in)
            for(l in 1:t){
              ll.y.right.cand[swapped,,l]=dbinom(y.right.tmp[swapped,,l],K[l],z[swapped,l]*(ones[swapped,,l]*pd1[swapped,,l]+twos[swapped,,l]*(2*pd1[swapped,,l]-pd1[swapped,,l]*pd1[swapped,,l])),log=TRUE)
            }
            if(!is.finite(sum(ll.y.right.cand))){
              stop("ll.y.right.cand not finite. Maybe swap.tol too large")
            }
            llswap.curr<- sum(ll.y.right[swapped,,l])
            llswap.cand<- sum(ll.y.right.cand[swapped,,l])
            lldiff=llswap.cand-llswap.curr
            #MH step
            if(runif(1)<exp(lldiff)*(jump.back/jump.probability) ){
              # lik.curr=lik.curr+lldiff #update likelihood
              y.right[swapped,,] <- y.right.tmp[swapped,,] #update right data
              ll.y.right[swapped,,]=ll.y.right.cand[swapped,,]
              ID_R<-newID #update data order
              map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
              zero.guys<- apply(y.both+y.left + y.right ,c(1,3),sum) == 0
              if(any(((colSums(y.both[guy1,,] + y.left[guy1,,]+ y.right.tmp[guy1,,])>0)==TRUE)&((z[guy1,]==0)==TRUE))){
                cat("error on guy1",fill=TRUE)
                browser()
              }
              if( any(((colSums(y.both[guy2,,] + y.left[guy2,,]+ y.right.tmp[guy2,,])>0)==TRUE)&((z[guy2,]==0)==TRUE))){
                cat("error on guy2",fill=TRUE)
                browser()
              }
            }
          }else{
            #joint update
            if(length(possible)>1){
              s.swap.in <-  sample( possible, 1)
            }
            if(length(possible) ==1){
              s.swap.in <- possible
            }
            if(length(possible)==0) next #no guys close enough to swap
            jump.probability<- 1/length(possible)  # h(theta*|theta)

            #compute  h(theta|theta*)
            trash<-   sqrt( (s1[s.swap.in,1]- s1[INcands,1])^2 + (s1[s.swap.in,2] - s1[INcands,2])^2  )
            trash<-  INcands[trash < swap.tol]
            jump.back<-  1/length(trash)

            ##  Which right encounter history is currently associated with both guy s.swap.in?
            guy2<- which(map[,2]==s.swap.in)
            if(guy1==guy2)next
            newID<-ID_R
            newID[guy1]<- s.swap.in
            newID[guy2]<- s.swap.out

            ## recompute 'data' and compute likelihood components for swapped guys only
            y.right.tmp<- apply( right[order(newID),3,,,], c(1,2,4), sum)#Only y.right changed
            swapped=c(s.swap.in,s.swap.out)

            ##update z before ll.y.right
            #Get likelihood for all possible z histories
            ll_z_possible[,1]=dbinom(zpossible[,1], 1, psi,log=TRUE)
            for(l in 2:t){
              Ezpossible[,l-1]=zpossible[,l-1]*phiuse[l-1] + apossible[,l-1]*gamma.prime[l-1]
              ll_z_possible[,l]=dbinom(zpossible[,l], 1, Ezpossible[,l-1],log=TRUE)
            }

            #Proposing y.right.tmp changes known.matrix
            swappedL=c(which(ID_L==swapped[1]),which(ID_L==swapped[2]))
            tmpdata.prop=both[swapped,,,,] + left[swappedL,,,,] + right[c(guy1,guy2),,,,]
            tmpdata.prop<- apply(tmpdata.prop,c(1,3,5),sum)
            known.matrix.prop=1*(apply(tmpdata.prop,c(1,3),sum)>0)

            #pick a new z based on updated known.matrix constraints
            zchoose=zcandguy=prop.prob=back.prob=rep(NA,2)
            zprop=matrix(NA,nrow=2,ncol=t)
            for(i2 in 1:2){
              #Zero out known matrix years for both swapped
              fixed=which(known.matrix.prop[i2,]==1)
              cancel=rep(1,nzpossible)
              for(i3 in 1:nzpossible){
                if(!all(fixed%in%which(zpossible[i3,]==1))){
                  cancel[i3]=0
                }
              }
              #new z stuff
              propto1=rowSums(exp(ll_z_possible))*(1*(cancel==1))
              propto=propto1/sum(propto1)
              zchoose[i2]=sample(1:nzpossible,1,prob=propto)
              zprop[i2,]=zpossible[zchoose[i2],]
              prop.prob[i2]=prod(propto[zchoose[i2]])
              #old z stuff
              zcandguy[i2]=which(apply(zpossible,1,function(x){all(x==z[swapped[i2],])}))
              back.prob[i2]=propto[zcandguy[i2]]
            }

            #Because a and z changes, must update gamma.prime and Ez
            #Don't need to update all years every time, but not figuring that out for now
            aprop=apossible[zchoose,]

            #Calculate gamma.prime, Ez, and ll.z candidates
            ll.z.cand[swapped,1] <- dbinom(zprop[,1], 1, psi, log=TRUE)
            z.cand=z
            z.cand[swapped,]=zprop
            a.cand=a
            a.cand[swapped,]=aprop
            Ntmp=colSums(z.cand)
            for(l in 2:t){
              gamma.prime.cand[l-1]=(Ntmp[l-1]*gammause[l-1]) / sum(a.cand[,l-1])
              if(gamma.prime.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
                warning("Rejected z due to low M")
                next
              }
              Ez.cand[,l-1]=z.cand[,l-1]*phiuse[l-1] + a.cand[,l-1]*gamma.prime.cand[l-1]
              ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
            }

            #update new ll.y.right for new y.right and z
            for(l in 1:t){
              ll.y.right.cand[swapped,,l]=dbinom(y.right.tmp[swapped,,l],K[l],zprop[,l]*(ones[swapped,,l]*pd1[swapped,,l]+twos[swapped,,l]*(2*pd1[swapped,,l]-pd1[swapped,,l]*pd1[swapped,,l])),log=TRUE)
            }
            if(!is.finite(sum(ll.y.right.cand))){
              stop("ll.y.right.cand not finite. Maybe swap.tol too large")
            }
            llyswap.curr<- sum(ll.y.right[swapped,,])
            llyswap.cand<- sum(ll.y.right.cand[swapped,,])
            llzswap.curr=sum(ll.z[swapped,1])+sum(ll.z[,2:t])
            llzswap.cand=sum(ll.z.cand[swapped,1])+sum(ll.z.cand[,2:t])
            #MH step
            if(runif(1)<exp(c(llyswap.cand+llzswap.cand)-c(llyswap.curr+llzswap.curr))*((jump.back*prod(back.prob))/(jump.probability*prod(prop.prob)))){
              y.right[swapped,,] <- y.right.tmp[swapped,,] #update right data
              ll.y.right[swapped,,]=ll.y.right.cand[swapped,,]
              ID_R<-newID #update data order
              map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
              known.matrix[swapped,]=known.matrix.prop
              gamma.prime=gamma.prime.cand
              N=Ntmp
              Ez=Ez.cand
              ll.z[swapped,1]=ll.z.cand[swapped,1]
              ll.z[,2:t]=ll.z.cand[,2:t]
              zero.guys<- apply(y.both+y.left + y.right ,c(1,3),sum) == 0
              if(any(((colSums(y.both[guy1,,] + y.left[guy1,,]+ y.right.tmp[guy1,,])>0)==TRUE)&((z[guy1,]==0)==TRUE))){
                cat("error on guy1",fill=TRUE)
                browser()
              }
              if( any(((colSums(y.both[guy2,,] + y.left[guy2,,]+ y.right.tmp[guy2,,])>0)==TRUE)&((z[guy2,]==0)==TRUE))){
                cat("error on guy2",fill=TRUE)
                browser()
              }
            }
          }
        }
        #Do any captured guys have z==0?
        if(any(!zero.guys&z==0)){
          cat("coming out of ID_R update: some z = 0!",fill=TRUE)
          browser()
        }
      }

      # #known matrix must be updated
      # tmpdata<- both + left[order(ID_L),,,,] + right[order(ID_R),,,,]
      # tmpdata<- apply(tmpdata,c(1,3,5),sum)
      # known.matrix=1*(apply(tmpdata,c(1,3),sum)>0)

      # Update z[,1]
      if(t>3){
        upz=which(!(z[,1]==0&z[,2]==0&rowSums(z[,3:t]>0))&known.matrix[,1]==0)
      }else if(t==3){
        upz=which(!(z[,1]==0&z[,2]==0&z[,3]>0)&known.matrix[,1]==0)
      }else{
        upz=which(known.matrix[,1]==0)
      }
      for(i in upz){
        gamma.prime.cand <- gamma.prime
        #Do we need to modify a in more than one year.
        z.cand <- z #use full z to calculate correct proposed Ez.cand
        z.cand[i,1] <- 1-z[i,1]
        ll.y.both.cand[i,,1]=dbinom(y.both[i,,1],K[1],z.cand[i,1]*pd2[i,,1],log=TRUE)
        ll.y.both.cand[i,which(X[[1]][,3]==1),1]=0 #cancel out contributions from single traps
        ll.y.left.cand[i,,1]=dbinom(y.left[i,,1],K[1],z.cand[i,1]*(ones[i,,1]*pd1[i,,1]+twos[i,,1]*(2*pd1[i,,1]-pd1[i,,1]*pd1[i,,1])),log=TRUE)
        ll.y.right.cand[i,,1]=dbinom(y.right[i,,1],K[1],z.cand[i,1]*(ones[i,,1]*pd1[i,,1]+twos[i,,1]*(2*pd1[i,,1]-pd1[i,,1]*pd1[i,,1])),log=TRUE)
        if(((z.cand[i,1]==1&sum(a[i,])==t)|(sum(z[i,])==1&sum(a[i,])==0))&(t>2)){#Are we turning on a guy that was never on before? or turning off a guy that was only on on z1?
          a.cand <- a
          #only a option is all on or all off
          if(z.cand[i,1]==1&sum(a[i,])==t){
            a.cand[i,]=0 #all off
          }else{
            a.cand[i,]=1 #all on
          }
          #Calculate gamma.prime, Ez, and ll.z candidates
          ll.z.cand[i,1] <- dbinom(z.cand[i,1], 1, psi, log=TRUE)
          #Make object to put N1_prop in to but use current values of the other Ns
          Ntmp=N
          Ntmp[1]=sum(z.cand[,1])
          for(l in 2:t){
            gamma.prime.cand[l-1]=(Ntmp[l-1]*gammause[l-1]) / sum(a.cand[,l-1])
            if(gamma.prime.cand[l-1] > 1) { # E(Recruits) must be < nAvailable
              warning("Rejected z due to low M")
              next
            }
            Ez.cand[,l-1]=z.cand[,l-1]*phiuse[l-1] + a.cand[,l-1]*gamma.prime.cand[l-1]
            ll.z.cand[,l]=dbinom(z.cand[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp((sum(ll.y.both.cand[i,,1])+sum(ll.y.left.cand[i,,1])+sum(ll.y.right.cand[i,,1])+ll.z.cand[i,1]+sum(ll.z.cand[,-1]))-(sum(ll.y.both[i,,1])+sum(ll.y.left[i,,1])+sum(ll.y.right[i,,1])+ll.z[i,1]+sum(ll.z[,-1])) )) {
            ll.y.both[i,,1] = ll.y.both.cand[i,,1]
            ll.y.left[i,,1] = ll.y.left.cand[i,,1]
            ll.y.right[i,,1] = ll.y.right.cand[i,,1]
            ll.z[i,1] = ll.z.cand[i,1]
            ll.z[,-1] = ll.z.cand[,-1]
            Ez = Ez.cand
            z[i,1] = z.cand[i,1]
            gamma.prime = gamma.prime.cand
            a=a.cand
          }
        }else{#Don't need to modify more than one year
          z1.cand <- z[,1]
          z1.cand[i] <- 1-z[i,1]
          a1.cand <- a[,1]
          a1.cand[i] <- 1-z1.cand[i]
          gamma.prime.cand[1] <- sum(z1.cand)*gamma[1] / sum(a1.cand)
          if(gamma.prime.cand[1] > 1) { # E(Recruits) must be < nAvailable
            warning("Rejected z due to low M")
            next
          }
          ll.z.cand[i,1] <- dbinom(z1.cand[i], 1, psi, log=TRUE)
          Ez.cand[,1]=z1.cand*phi[1] + a1.cand*gamma.prime.cand[1]
          ll.z.cand[,2] <- dbinom(z[,2], 1, Ez.cand[,1], log=TRUE)
          if(runif(1) < exp((sum(ll.y.both.cand[i,,1])+sum(ll.y.left.cand[i,,1])+sum(ll.y.right.cand[i,,1])+ ll.z.cand[i,1]+sum(ll.z.cand[,2]))-(sum(ll.y.both[i,,1])+sum(ll.y.left[i,,1])+sum(ll.y.right[i,,1])+ll.z[i,1]+sum(ll.z[,2])) )) {#z1 and z2 matter
            ll.y.both[i,,1] = ll.y.both.cand[i,,1]
            ll.y.left[i,,1] = ll.y.left.cand[i,,1]
            ll.y.right[i,,1] = ll.y.right.cand[i,,1]
            ll.z[i,1] = ll.z.cand[i,1]
            ll.z[,2] = ll.z.cand[,2]
            Ez[,1] = Ez.cand[,1]
            z[i,1] = z1.cand[i]
            gamma.prime[1] = gamma.prime.cand[1]
            a[i,1]=a1.cand[i]
          }
        }
      }
      N[1]=sum(z[,1])

      for(l in 2:t){
        #Determine who can be updated
        upz=which(known.matrix[,l]==0)#guys not caught on or on either side of this occ
        #Figure out illegal moves.  Depends on t.
        #Could use apply function, but need something to convert to Rcpp
        if(t>3){
          if(l==2&l==(t-2)){#only happens if t=4
            rem=which(z[,1]>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
            rem2=which(z[,l]==0&z[,l+1]==0&z[,t]>0)#guys with 0 0 and subsequent 1 can't be turned on
          }else if(l==2){
            rem=which(z[,1]>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
            rem2=which(z[,l]==0&z[,l+1]==0&rowSums(z[,(l+2):t])>0) #guys with 0 0 and subsequent 1 can't be turned on
          }else if(l==(t-2)){
            rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
            rem2=which(z[,l]==0&z[,l+1]==0&z[,t]>0)#guys with 0 0 and subsequent 1 can't be turned on
          }else if(l==(t-1)){
            rem=which(rowSums(z[,1:(l-1)])>0&z[,t]>0)#guys in pop before and after
            rem2=integer()
          }else if (l==t){
            rem=integer()
            rem2=integer()
          }else{
            rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)#guys in pop before and after
            rem2=which(z[,l]==0&z[,l+1]==0&rowSums(z[,(l+2):t])>0)#guys with 0 0 and subsequent 1 can't be turned on
          }
        }else if(t==3){
          if(l==2){
            rem=which(z[,1]>0&z[,t]>0)#guys in pop before and after
            rem2=integer()
          }else{#l==3
            rem=integer()
            rem2=integer()
          }
        }else{#t==2
          rem=integer()
          rem2=integer()
        }
        rem3=which( (z[,l-1]==0) & (a[,l-1]==0)) #dead guys.
        remall=c(rem,rem2,rem3)
        upz=upz[!upz%in%remall]#We can change these guys
        navail=length(upz)
        if(navail<1){
          next
        }
        if(navail < proppars$propz[l-1]) {
          propz=navail
          warning("M isn't big enough to propose all propz")
        }else{
          propz=proppars$propz[l-1]
        }
        swapz=upz[sample.int(navail, propz)]
        #Update swapz one at a time
        for(i in swapz){
          zt.cand=z[,l]
          ####Try proposing not based on Ez. Just swap it. Take out prop and back probs. Should help mixing
          #Normal stuff
          zt.cand[i]=1-z[i,l]
          at.cand=1*(a[,l-1]==1&zt.cand==0) #who was available on last occasion and not proposed to be captured?
          ll.y.left.cand[i,,l]=dbinom(y.left[i,,l],K[l],zt.cand[i]*(ones[i,,l]*pd1[i,,l]+twos[i,,l]*(2*pd1[i,,l]-pd1[i,,l]*pd1[i,,l])),log=TRUE)
          ll.y.right.cand[i,,l]=dbinom(y.right[i,,l],K[l],zt.cand[i]*(ones[i,,l]*pd1[i,,l]+twos[i,,l]*(2*pd1[i,,l]-pd1[i,,l]*pd1[i,,l])),log=TRUE)
          ll.y.both.cand[i,,l]=dbinom(y.both[i,,l],K[l],zt.cand[i]*pd2[i,,l],log=TRUE)
          ll.y.both.cand[i,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
          # ll.z.cand[,l] <- dbinom(zt.cand, 1, Ez[,l-1], log=TRUE) ## Don't subset z
          ll.z.cand[i,l] <- dbinom(zt.cand[i], 1, Ez[i,l-1], log=TRUE) ## why not?
          # prior.z <- sum(ll.z[i,l])
          # prior.z.cand <- sum(ll.z.cand[i,l])
          prior.z=ll.z[i,l]
          prior.z.cand=ll.z.cand[i,l]
          #New stuff
          fix1=zt.cand[i]==1&sum(z[i,])==0 #guys never in pop proposed to be turned on
          fix2=sum(z[i,])==1&zt.cand[i]==0&z[i,l]==1 #guys in pop only once and proposed to be turned off
          if((fix1|fix2)&(t>3)&(l<t)){
            a.cand <- a
            a.cand[,l]=at.cand
            z.cand=z
            z.cand[,l]=zt.cand
            if(fix1){
              a.cand[i,l:t]=0 #all l:t off
            }
            if(fix2){
              a.cand[i,l:t]=1 #all on 1:(l-1) are already on
            }
            #Calculate gamma.prime, Ez, and ll.z for l:t
            reject=FALSE
            Ntmp=N
            Ntmp[l]=sum(zt.cand)
            for(l2 in l:(t-1)){
              gamma.prime.cand[l2]=(Ntmp[l2]*gammause[l2]) / sum(a.cand[,l2])
              if(gamma.prime.cand[l2] > 1){
                reject=TRUE
              }
              Ez.cand[,l2]=z.cand[,l2]*phiuse[l2] + a.cand[,l2]*gamma.prime.cand[l2]
              ll.z.cand[,l2+1]=dbinom(z.cand[,l2+1], 1, Ez.cand[,l2], log=TRUE)
              prior.z <- prior.z + sum(ll.z[,l2+1])
              prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l2+1])
            }
            if(reject){
              warning("Rejected z due to low M")
              next
            }
          }else{
            #Calculate gamma.prime, Ez, and ll.z for l
            if(l<t){ ## NOTE: Don't subset with swapz
              gamma.prime.cand[l] <- sum(zt.cand)*gammause[l] / sum(at.cand)
              if(gamma.prime.cand[l] > 1){
                warning("Rejected z due to low M")
                next
              }
              Ez.cand[,l] <- zt.cand*phiuse[l] + at.cand*gamma.prime.cand[l]
              ll.z.cand[,l+1] <- dbinom(z[,l+1], 1, Ez.cand[,l], log=TRUE)
              prior.z <- prior.z + sum(ll.z[,l+1])
              prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l+1])
            }
          }
          if(runif(1) < exp((sum(ll.y.both.cand[i,,l])+sum(ll.y.left.cand[i,,l])+sum(ll.y.right.cand[i,,l]) + prior.z.cand) - (sum(ll.y.both[i,,l])+sum(ll.y.left[i,,l])+sum(ll.y.right[i,,l]) +prior.z) )) {
            ll.y.both[i,,l] <- ll.y.both.cand[i,,l]
            ll.y.left[i,,l] <- ll.y.left.cand[i,,l]
            ll.y.right[i,,l] <- ll.y.right.cand[i,,l]
            # ll.z[,l] <- ll.z.cand[,l]
            ll.z[i,l] <- ll.z.cand[i,l]
            if((fix1|fix2)&(t>3)&(l<t)){
              z=z.cand
              a=a.cand
              ll.z[,l:t]=ll.z.cand[,l:t]
              Ez[,l:(t-1)]=Ez.cand[,l:(t-1)]
              gamma.prime[l:(t-1)]= gamma.prime.cand[l:(t-1)]
            }else{
              z[,l] <- zt.cand
              a[,l] <- at.cand
              if(l < t) {
                ll.z[,l+1] <- ll.z.cand[,l+1]
                Ez[,l] <- Ez.cand[,l]
                gamma.prime[l] <- gamma.prime.cand[l]
              }
            }
            N[l] <- sum(z[,l])
          }
        }
      }
      if(t>3){
        #a for last year isn't updated so it does not matter if it comes back on.
        sanity=any(rowSums(z)==0&rowSums(a[,-t])!=(t-1))
        for(l2 in 2:(t-2)){
          sanity=c(sanity,any(a[,l2]==0&a[,l2-1]==1&a[,l2+1]==1))
        }
        if(any(sanity)){stop("insanity")}
      }
      #update psi
      psi <- rbeta(1, 1+N[1], 1+M-N[1])
      ll.z[,1] <- ll.z.cand[,1] <- dbinom(z[,1], 1, psi, log=TRUE)

      #Update phi
      if(length(phi)==(t-1)){#if time-specific survival
        for(l in 2:t){
          survive=sum(z[,l-1]==1&z[,l]==1)
          dead=sum(z[,l-1]==1&z[,l]==0)
          phi[l-1]=rbeta(1, 1+survive, 1+dead)
        }
      }else{
        survive=sum(z[,-t]==1&z[,-1]==1)
        dead=sum(z[,-t]==1&z[,-1]==0)
        phi=rbeta(1, 1+survive, 1+dead)
      }

      ## Update gamma
      ## NOTE: Must update ll.z, Ez, etc...
      if(length(phi)==1){
        phiuse=rep(phi,t-1)
      }else{
        phiuse=phi
      }
      if(length(gamma)==1){
        gamma.cand <- rnorm(1, gamma, proppars$gamma)
        gamma.cand.ok <- TRUE
        for(l in 2:t) {
          gamma.prime.cand[l-1] <- (N[l-1]*gamma.cand) / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=FALSE
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases. originally. Changed from Chandler. changed back.
          for(l in 2:t) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            # ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
          }
          if(runif(1) < exp(sum(ll.z.cand[,-1]) - sum(ll.z[,-1]))) {
            gamma <- gamma.cand
            gamma.prime <- gamma.prime.cand
            Ez <- Ez.cand
            ll.z[,-1] <- ll.z.cand[,-1]
          }
        }
      }else{
        for(l in 2:t){
          gamma.cand <- rnorm(1, gamma[l-1], proppars$gamma[l-1])
          gamma.cand.ok <- TRUE
          gamma.prime.cand[l-1] <- N[l-1]*gamma.cand / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=!gamma.cand.ok
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
          if(gamma.cand>0 & gamma.cand.ok) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            # ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
            ll.z.cand[,l] <- dbinom(z[,l], 1, Ez.cand[,l-1], log=TRUE)
            if(runif(1) < exp(sum(ll.z.cand[,l]) - sum(ll.z[,l]))) {
              gamma[l-1] <- gamma.cand
              gamma.prime[l-1] <- gamma.prime.cand[l-1]
              Ez[,l-1] <- Ez.cand[,l-1]
              ll.z[,l] <- ll.z.cand[,l]
            }
          }
        }
      }
      #Update gamma use
      if(length(gamma)==1){
        gammause=rep(gamma,t-1)
      }else{
        gammause=gamma
      }
      ## Now we have to update the activity centers
      if(metamu){
        #Update within year ACs
        for (i in 1:M){
          for(l in 1:t){
            Scand=c(rnorm(1, s2[i,l,1], proppars$s2x), rnorm(1, s2[i,l,2], proppars$s2y))
            if(useverts==FALSE){
              inbox=Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
            }else{
              inbox=inout(Scand,vertices)
            }
            if(inbox) {
              dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
              if(length(lam01)==1&length(sigma==1)){
                lamd1.cand[i,1:nrow(X[[l]]),l]<- lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
                lamd2.cand[i,1:nrow(X[[l]]),l]<- lam02*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else if(length(lam01)==t&length(sigma)==1){
                lamd1.cand[i,1:nrow(X[[l]]),l]<- lam01[l]*exp(-dtmp*dtmp/(2*sigma*sigma))
                lamd2.cand[i,1:nrow(X[[l]]),l]<- lam02[l]*exp(-dtmp*dtmp/(2*sigma*sigma))
              }else{
                lamd1.cand[i,1:nrow(X[[l]]),l]<- lam01[l]*exp(-dtmp*dtmp/(2*sigma[l]*sigma[l]))
                lamd2.cand[i,1:nrow(X[[l]]),l]<- lam02[l]*exp(-dtmp*dtmp/(2*sigma[l]*sigma[l]))
              }
              pd1.cand[i,,l]=1-exp(-lamd1.cand[i,,l])
              pd2.cand[i,,l]=1-exp(-lamd2.cand[i,,l])
              ll.y.left.cand[i,,l]=dbinom(y.left[i,,l],K[l],z[i,l]*(ones[i,,l]*pd1.cand[i,,l]+twos[i,,l]*(2*pd1.cand[i,,l]-pd1.cand[i,,l]*pd1.cand[i,,l])),log=TRUE)
              ll.y.right.cand[i,,l]=dbinom(y.right[i,,l],K[l],z[i,l]*(ones[i,,l]*pd1.cand[i,,l]+twos[i,,l]*(2*pd1.cand[i,,l]-pd1.cand[i,,l]*pd1.cand[i,,l])),log=TRUE)
              ll.y.both.cand[i,,l]=dbinom(y.both[i,,l],K[l],pd2.cand[i,,l]*z[i,l],log=TRUE)
              ll.y.both.cand[i,which(X[[l]][,3]==1),l]=0 #cancel out contributions from single traps
              ll.s2.cand<- dnorm(Scand[1],s1[i,1],sigma_t,log=TRUE)+dnorm(Scand[2],s1[i,2],sigma_t,log=TRUE)
              if(runif(1) < exp((sum(ll.y.both.cand[i,,l])+sum(ll.y.left.cand[i,,l])+sum(ll.y.right.cand[i,,l])+ll.s2.cand) -(sum(ll.y.both[i,,l])+sum(ll.y.left[i,,l])+sum(ll.y.right[i,,l])+ll.s2[i,l]))){
                s2[i,l,] <- Scand
                D[i,,l] <- dtmp
                lamd1[i,,l] <- lamd1.cand[i,,l]
                lamd2[i,,l] <- lamd2.cand[i,,l]
                pd1[i,,l]=pd1.cand[i,,l]
                pd2[i,,l]=pd2.cand[i,,l]
                ll.y.both[i,,l]=ll.y.both.cand[i,,l]
                ll.y.left[i,,l]=ll.y.left.cand[i,,l]
                ll.y.right[i,,l]=ll.y.right.cand[i,,l]
                ll.s2[i,l]=ll.s2.cand
              }
            }
          }
        }
        #Update meta mus
        for (i in 1:M){
          Scand <- c(rnorm(1,s1[i,1],proppars$s1x), rnorm(1,s1[i,2],proppars$s1y))
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            inbox=inout(Scand,vertices)
          }
          if(inbox){
            #Count z==0 guys
            ll.s2.cand<-dnorm(s2[i,,1],Scand[1],sigma_t,log=TRUE)+dnorm(s2[i,,2],Scand[2],sigma_t,log=TRUE)
            if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2[i,]))) {
              s1[i, ]=Scand
              ll.s2[i,]=ll.s2.cand
            }
          }
        }
        #Update sigma_t
        sigma_t.cand <- rnorm(1,sigma_t,proppars$sigma_t)
        if(sigma_t.cand > 0){
          ll.s2.cand=dnorm(s2[,,1],s1[,1],sigma_t.cand,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t.cand,log=TRUE)
          if (runif(1) < exp(sum(ll.s2.cand) - sum(ll.s2))) {
            sigma_t=sigma_t.cand
            ll.s2=ll.s2.cand
          }
        }
      }else{#Stationary ACs
        for (i in 1:M) {
          Scand <- c(rnorm(1, s1[i, 1], proppars$s2x), rnorm(1, s1[i, 2], proppars$s2y))
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            inbox=inout(Scand,vertices)
          }
          if (inbox) {
            dtmp=matrix(Inf,maxJ,t)
            for(l in 1:t){
              dtmp[1:nrow(X[[l]]),l] <- sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
            }
            if(length(lam01)==1&length(sigma==1)){
              lamd1.cand[i,,]<- lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
              lamd2.cand[i,,]<- lam02*exp(-dtmp*dtmp/(2*sigma*sigma))
            }else if(length(lam01)==t&length(sigma)==1){
              for(l in 1:t){
                lamd1.cand[i,,l]<- lam01[l]*exp(-dtmp[,l]*dtmp[,l]/(2*sigma*sigma))
                lamd2.cand[i,,l]<- lam02[l]*exp(-dtmp[,l]*dtmp[,l]/(2*sigma*sigma))
              }
            }else{
              for(l in 1:t){
                lamd1.cand[i,,l]<- lam01[l]*exp(-dtmp[,l]*dtmp[,l]/(2*sigma[l]*sigma[l]))
                lamd2.cand[i,,l]<- lam02[l]*exp(-dtmp[,l]*dtmp[,l]/(2*sigma[l]*sigma[l]))
              }
            }
            pd1.cand[i,,]=1-exp(-lamd1.cand[i,,])
            pd2.cand[i,,]=1-exp(-lamd2.cand[i,,])
            ll.y.both.cand[i,,]=ll.y.both[i,,]
            ll.y.left.cand[i,,]=ll.y.left[i,,]
            ll.y.right.cand[i,,]=ll.y.right[i,,]
            for(l in 1:t) {
              if(z[i,l]==0)
                next
              ll.y.left.cand[i,,l]=dbinom(y.left[i,,l],K[l],z[i,l]*(ones[i,,l]*pd1.cand[i,,l]+twos[i,,l]*(2*pd1.cand[i,,l]-pd1.cand[i,,l]*pd1.cand[i,,l])),log=TRUE)
              ll.y.right.cand[i,,l]=dbinom(y.right[i,,l],K[l],z[i,l]*(ones[i,,l]*pd1.cand[i,,l]+twos[i,,l]*(2*pd1.cand[i,,l]-pd1.cand[i,,l]*pd1.cand[i,,l])),log=TRUE)
              ll.y.both.cand[i,,l]=dbinom(y.both[i,,l],K[l],z[i,l]*pd2.cand[i,,l],log=TRUE)
              ll.y.both.cand[i,which(X[[l]][,3]==1),l]=0
            }
            if(runif(1) < exp((sum(ll.y.both.cand[i,,])+sum(ll.y.left.cand[i,,])+sum(ll.y.right.cand[i,,])) -(sum(ll.y.both[i,,])+sum(ll.y.left[i,,])+sum(ll.y.right[i,,])))){
              s1[i, ] <- Scand
              D[i,, ] <- dtmp
              lamd1[i,, ] <- lamd1.cand[i,,]
              lamd2[i,, ] <- lamd2.cand[i,,]
              pd1[i,,]=pd1.cand[i,,]
              pd2[i,,]=pd2.cand[i,,]
              ll.y.both[i,,]=ll.y.both.cand[i,,]
              ll.y.left[i,,]=ll.y.left.cand[i,,]
              ll.y.right[i,,]=ll.y.right.cand[i,,]
            }
          }
        }
      }

      #Do we record output on this iteration?
      if(iter>(nburn)&iter%%nthin==0){
        s1xout[idx,]<- s1[,1]
        s1yout[idx,]<- s1[,2]
        zout[idx,,]<- z
        ID_Lout[idx,]=ID_L
        ID_Rout[idx,]=ID_R
        if(metamu){
          out[idx,]<- c(lam01,lam02,sigma ,gamma,phi,N,sigma_t)
          s2xout[idx,,]<- s2[,,1]
          s2yout[idx,,]<- s2[,,2]
        }else{
          out[idx,]<- c(lam01,lam02,sigma ,gamma,phi,N)
        }
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      if(metamu){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout,ID_Lout=ID_Lout,ID_Rout=ID_Rout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout,ID_Lout=ID_Lout,ID_Rout=ID_Rout)
      }
    }else{
      list(out=out)
    }
  }

