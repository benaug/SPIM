SCRmcmcOpenSPIMRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,swap=10,swap.tol=1,
           proppars=list(lam01=0.05,lam01=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,jointZ=TRUE){
    library(abind)
    t=dim(data$both)[4]
    both<-data$both
    left<-data$left
    right<-data$right
    if(dim(right)[1]>dim(left)[1]){
      storeL=left
      storeR=right
      dimL=dim(left)
      left=array(0,dim=dim(right))
      right=array(0,dim=dimL)
      left=storeR
      right=storeL
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
      stop("Input either 1 or t initial values for lam01")
    }
    if(!length(lam02)%in%c(1,t)){
      stop("Input either 1 or t initial values for lam02")
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
    if((jointZ==FALSE)&length(proppars$propz)!=(t-1)){
      stop("must supply t-1 proppars for propz")
    }
    if(length(lam01)!=length(proppars$lam01)){
      stop("Must supply a tuning parameter for each lam01")
    }
    if(length(lam02)!=length(proppars$lam02)){
      stop("Must supply a tuning parameter for each lam02")
    }
    if(length(sigma)!=length(proppars$sigma)){
      stop("Must supply a tuning parameter for each sigma")
    }
    if(length(gamma)!=length(proppars$gamma)){
      stop("Must supply a tuning parameter for each gamma")
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
    both<- abind(both,array(0, dim=c( M-dim(both)[1],maxJ,maxK, t)), along=1)
    left<- abind(left,array(0, dim=c( M-dim(left)[1],maxJ,maxK, t)), along=1)
    right<- abind(right,array(0, dim=c( M-dim(right)[1],maxJ,maxK, t)), along=1)

    #Make initial complete data set
    tmpdata<- both + left[order(ID_L),,,] + right[order(ID_R),,,]
    tmpdata<- apply(tmpdata,c(1,2,4),sum)
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
      nremain=nrecruit-length(recruits)
      if(nremain>0){
        cands=setdiff(which(a[,i-1]==1),recruits)#Who is available to recruit?
        cands=cands[!cands%in%knownguys]#Don't mess with known guys
        if(nremain>length(cands)){
          nremain=length(cands)
        }
        pick=sample(cands,nremain)
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Can't initialize all E[recruits] given initial value for gamma. Should probably raise M?")
      }
      #Survival
      nlive=rbinom(1,sum(z[,i-1]),phisim[i-1])#How many to live
      livers=which(z[,i-1]==1&z[,i]==1)
      nremain=nlive-length(livers)
      a[livers,i]=0
      if(nremain>0){
        cands=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1)
        cands=cands[!cands%in%knownguys]
        if(nremain>length(cands)){
          nremain=length(cands)
        }
        pick=sample(cands,nremain)
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
    }else{
      proppars$s1x=proppars$s1y=0
      ll.s2=matrix(NA,nrow=M,ncol=t)
      sigma_t=0
      proppars$sigma_t=0
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
      colnames(out)<-c(lam01names,lam02names,sigmanames,gammanames,phinames,Nnames,"sigma_t")
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

    y.both<- apply(both[,,,], c(1,2,4), sum)
    y.left<- apply(left[order(ID_L),,,], c(1,2,4), sum)
    y.right<- apply(right[order(ID_R),,,], c(1,2,4), sum)
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
    if(!is.finite(ll.y.left.sum)|!is.finite(ll.y.right.sum)){
      stop("Single side detection function starting values produce -Inf log likelihood values. Try increasing sigma and/or lam01")
    }

    #Figure out all possible z histories
    zpossible=cbind(c(1,1,0),c(1,0,1))
    if (t > 2) {#4 lines from From Rcapture histpost()
      for (i in (3:t)) {
        zpossible=cbind(c(rep(1,2^(i-1)),rep(0,((2^(i-1))-1))), rbind(zpossible, rep(0, (i-1)), zpossible))
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
    }

    zpossible=rbind(zpossible,rep(0,t))#add on all zero history
    nzpossible=nrow(zpossible)
    apossible=matrix(1,nrow=nzpossible,ncol=t)
    apossible[zpossible[,1]==1,1]=0
    for(l in 2:t){
      apossible[apossible[,l-1]==1&zpossible[,l]==1,l]=0
      apossible[apossible[,l-1]==0,l]=0
    }
    Ezpossible=matrix(NA,nrow=nzpossible,ncol=t-1)
    ll_z_possible=matrix(NA,nrow=nzpossible,ncol=t)
    #Can always update anyone who wasn't captured on every occasion
    upz3=which(rowSums(known.matrix)!=t)
    #Zero out known matrix years for both swapped
    cancel=matrix(1,nrow=M,ncol=nzpossible)
    for(i in 1:M){
      fixed=which(known.matrix[i,]==1)
      for(i2 in 1:nzpossible){
        if(!all(fixed%in%which(zpossible[i2,]==1))){
          cancel[i,i2]=0
        }
      }
    }

    Xidx=unlist(lapply(X,dim))[seq(1,2*length(X)-1,2)]
    Xcpp=array(NA,dim=c(t,max(Xidx),3))
    for(l in 1:t){
      for(j in 1:Xidx[l]){
        Xcpp[l,j,]=as.numeric(X[[l]][j,])
      }
    }
    each=unlist(lapply(inits,length))[1:5]
    npar=sum(each)+t
    if(metamu){
      npar=npar+1
    }
    ones=1*ones[1,,]
    twos=1*twos[1,,]
    #So these aren't modified by Rcpp
    lam01in=lam01
    lam02in=lam02
    sigmain=sigma
    gammain=gamma
    phiin=phi
    updates=c(uplam01,uplam02,nleft>0,nright>0)
    #only need 3d left and right if no tf
    # both=apply(both,c(1,2,4),sum)
    left=apply(left,c(1,2,4),sum)
    right=apply(right,c(1,2,4),sum)

    store=mcmc_Open_SPIM(  lam01in,lam02in,  sigmain,  gammain, gamma.prime,  phiin, D,lamd1,lamd2, y.both,y.left,y.right, z, a,  s1, s2,
                           metamu, useverts, vertices, xlim, ylim, known.matrix, Xidx, Xcpp, K, Ez,  psi,N, proppars$lam01,proppars$lam02,
                           proppars$sigma, proppars$propz,  proppars$gamma, proppars$s1x,  proppars$s1y,proppars$s2x,proppars$s2y,
                           proppars$sigma_t,sigma_t,niter,nburn,nthin,npar,each,jointZ, zpossible,apossible,cancel,ID_L,ID_R,ones,twos,
                           updates,swap,swap.tol,Nfixed,left,right)
    out=store[[1]]
    s1xout=store[[2]]
    s1yout=store[[3]]
    s2xout=store[[4]]
    s2yout=store[[5]]
    zout=store[[6]]
    ID_Lout=store[[7]]
    ID_Rout=store[[8]]
    warn=store[[9]]


    known.matrix=store[[17]]
    cancel=matrix(1,nrow=M,ncol=nzpossible)
    for(i in 1:M){
      fixed=which(known.matrix[i,]==1)
      for(i2 in 1:nzpossible){
        if(!all(fixed%in%which(zpossible[i2,]==1))){
          cancel[i,i2]=0
        }
      }
    }
    cancel==store[[16]]


    # llyLcandsum=store[[10]]
    # llyLsum=store[[11]]
    # llzcandsum=store[[12]]
    # llzsum=store[[13]]
    # backprobZID=store[[14]]
    # propprobZID=store[[15]]
    # backprob=store[[16]]
    # propprob=store[[17]]

    # exp((llyLcandsum+llzcandsum)-(llyLsum+llzsum))*((backprob*prod(backprobZID))/(propprob*prod(propprobZID)))
    # store[[18]]
    # store[[19]]
    # store[[20]]
    # store[[21]]
    # knownmatrixprop=store[[10]]
    # swapped=store[[11]]
    # swappedR=store[[12]]
    # guys=store[[13]]
    # cancelp=store[[14]]
    #
    # cancel.prop=matrix(1,nrow=2,ncol=nzpossible)
    # for(i2 in 1:2){
    #   #Zero out known matrix years for both swapped
    #   fixed=which(knownmatrixprop[i2,]==1)
    #   for(i3 in 1:nzpossible){
    #     if(!all(fixed%in%which(zpossible[i3,]==1))){
    #       cancel.prop[i2,i3]=0
    #     }
    #   }
    # }
    # cancelp
    # cancel.prop




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


