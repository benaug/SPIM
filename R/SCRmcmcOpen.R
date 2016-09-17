SCRmcmcOpen <-
  function(data,niter=2400,nburn=1200, nthin=5,M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    library(abind)
    y<-data$y
    X<-data$X
    J<-data$J
    maxJ=max(J)
    K<-data$K
    n=data$n
    n2d<-data$n2d
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
    # psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    sigma_t=inits$sigma_t
    gamma=inits$gamma
    phi=inits$phi
    psi=inits$psi
    if(!length(lam0)%in%c(1,t)){
      stop("Input either 1 or t initial values for lam0")
    }
    if(!length(sigma)%in%c(1,t)){
      stop("Input either 1 or t initial values for sigma")
    }
    if(!length(gamma)%in%c(1,t-1)){
      stop("Input either 2 or t initial values for gammma (psi and fixed gamma)")
    }
    if(!length(phi)%in%c(1,t-1)){
      stop("Input either 1 or t initial values for phi")
    }

    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],maxJ, t)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))

    #Initialize z, r, and a consistent with y
    known.matrix=1*(apply(y,c(1,3),sum)>0)#make z consistent with y
    for(l in 2:(t-1)){#Turn on zeros with 1's on either side
      known.matrix[known.matrix[,l]==0&known.matrix[,l-1]==1&rowSums(matrix(known.matrix[,(l+1):t],nrow=M)>0),l]=1
    }
    z=known.matrix
    # r=array(0,dim=dim(z))
    #turn on z's with caps on either side so we know they were in pop
    for(i in 1:n){
      idx=which(z[i,]==1)
      z[i,min(idx):max(idx)]=1
    }
    z[(n+1):M,1]=rbinom(M-n,1,psi)#add augmented guys to t=1 with psi
    # r[,1]=z[,1]
    a=matrix(1,nrow=M,ncol=t) #a is available to be recruited
    # a[,1]=1-z[,1]
    a[which(z[,1]==1),]=0
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
      nrecruit=round(sum(z[,i-1])*gammasim[i-1])
      recruits=which(z[1:n,i-1]==0&z[1:n,i]==1)#who is already recruited based on y constraints?
      a[which(z[,i]==1),i:t]=0 #anyone alive in time t can't be recruited
      nleft=nrecruit-length(recruits)
      if(nleft>0){
        cands=setdiff(which(a[,i-1]==1),recruits)#Who is available to recruit?
        cands=cands[!cands%in%1:n]#Don't mess with known guys
        if(nleft>length(cands)){
          nleft=length(cands)
        }
        pick=sample(cands,nleft)
        z[pick,i]=1
        a[pick,i:t]=0 #no longer available for recruit on any occasion
      }else{
        warning("Should probably raise M")
      }
      #Survival
      nlive=rbinom(1,sum(z[,i-1]),phisim[i-1])#How many to live
      livers=which(z[1:n,i-1]==1&z[1:n,i]==1)
      nleft=nlive-length(livers)
      a[livers,i]=0
      if(nleft>0){
        cands=which(apply(matrix(z[,i:t],nrow=M),1,sum)==0&z[,i-1]==1)
        cands=cands[!cands%in%1:n]
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
    ll.z.cand=ll.z
    #Optimize starting locations given where they are trapped. Initalizing s1 and s2 at same locs
    s1<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(known.vector==1) #switch for those actually caught
    for(i in idx){
      trps=matrix(0,nrow=0,ncol=2)
      for(j in 1:t){ #loop over t to get all cap locs
        trps<- rbind(trps,X[[j]][y[i,,j]>0,1:2])
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
        idx=which(rowSums(y[,,l])>0) #switch for those actually caught
        for(i in 1:M){
          if(i%in%idx){
            trps<- X[[l]][y[i,,l]>0,1:2]
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
      ll.s2=dnorm(s2[,,1],s1[,1],sigma_t,log=TRUE)+dnorm(s2[,,2],s1[,2],sigma_t,log=TRUE)
    }
    # some objects to hold the MCMC simulation output
    if(niter<(nburn)){
      stop("niter is smaller than nburn")
    }
    nstore=(niter-nburn)/nthin
    if((nburn)%%nthin!=0){
      nstore=nstore+1
    }
    if(length(lam0)==t){
      lam0names=paste("lam0",1:t,sep="")
    }else{
      lam0names="lam0"
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
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t+1)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames,"sigma_t")
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
      s2xout<- s2yout<-array(NA,dim=c(nstore,M,t))
    }else{
      out<-matrix(NA,nrow=nstore,ncol=length(lam0)+length(sigma)+length(gamma)+length(phi)+t)
      colnames(out)<-c(lam0names,sigmanames,gammanames,phinames,Nnames)
      s1xout<- s1yout<- matrix(NA,nrow=nstore,ncol=M)
      zout<-array(NA,dim=c(nstore,M,t))
    }
    idx=1 #for storing output not recorded every iteration

    D=lamd=ll.y=ll.y.cand=array(NA,dim=c(M,maxJ,t))
    D[is.na(D)]=Inf  #hack to allow years with different J and K to fit in one array
    for(i in 1:t){
      D[,1:nrow(X[[i]]),i]=e2dist(s2[,i,],X[[i]])
      if(length(lam0)==t&length(sigma)==t){
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam0)==1&length(sigma)==t){
        lamd[,,i]=lam0*exp(-D[,,i]^2/(2*sigma[i]*sigma[i]))
      }else if(length(lam0)==t&length(sigma)==1){
        lamd[,,i]=lam0[i]*exp(-D[,,i]^2/(2*sigma*sigma))
      }else{
        lamd[,,i]=lam0*exp(-D[,,i]^2/(2*sigma*sigma))
      }
    }
    #Calculate ll for observation model
    pd=pd.cand=1-exp(-lamd)
    for(l in 1:t){
      ll.y[,,l]= dbinom(y[,,l],K[l],pd[,,l]*z[,l],log=TRUE)
    }
    ll.y.cand=ll.y
    ll.y.t.sum=ll.y.cand.t.sum=apply(ll.y,3,sum) #ll summed for each year
    ll.y.sum=sum(ll.y.t.sum) #full ll sum
    lamd.cand=lamd
    #Check proppars
    if(length(proppars$propz)!=(t-1)){
      stop("must supply t-1 proppars for propz")
    }


    for(iter in 1:niter){
      # Update lam0
      if(length(lam0)==t){ #if lam0 is year-specific
        ll.y.t.sum=apply(ll.y,3,sum) #only needed for detection parameters, changes in z and AC updates
        for(l in 1:t){
          lam0.cand<- rnorm(1,lam0[l],proppars$lam0)
          if(lam0.cand > 0){
            if(length(sigma)==t){#if sigma is year specific
              lamd.cand[,,l]<- lam0.cand*exp(-D[,,l]^2/(2*sigma[l]*sigma[l]))
            }else{#fixed sigma
              lamd.cand[,,l]<- lam0.cand*exp(-D[,,l]^2/(2*sigma*sigma))
            }
            pd.cand[,,l]=1-exp(-lamd.cand[,,l])
            ll.y.cand[,,l]= dbinom(y[,,l],K[l],pd.cand[,,l]*z[,l],log=TRUE) #only need to update this year
            ll.y.cand.t.sum[l]=sum(ll.y.cand[,,l])#just 1 year
            if(runif(1) < exp(ll.y.cand.t.sum[l] -ll.y.t.sum[l])){
              lam0[l]<- lam0.cand
              lamd[,,l]=lamd.cand[,,l]
              pd[,,l]=pd.cand[,,l]
              ll.y[,,l]=ll.y.cand[,,l]
              ll.y.t.sum[l]=ll.y.cand.t.sum[l]
            }
            ll.y.sum=sum(ll.y.t.sum)
          }
        }
      }else{#fixed lam0
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          ll.y.sum=sum(ll.y)#Don't think I need to recalculate
          lamd.cand<- lam0.cand*exp(-D^2/(2*sigma*sigma)) #works for either single or multiple sigmas
          pd.cand=1-exp(-lamd.cand)
          for(l in 1:t){
            ll.y.cand[,,l]= dbinom(y[,,l],K[l],pd.cand[,,l]*z[,l],log=TRUE)
          }
          ll.y.cand.sum=sum(ll.y.cand)
          if(runif(1) < exp(ll.y.cand.sum - ll.y.sum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
            ll.y.sum=ll.y.cand.sum
          }
        }
      }
      #Update sigma
      if(length(sigma)==t){ #if sigma is year-specific
        for(l in 1:t){
          sigma.cand<- rnorm(1,sigma[l],proppars$sigma)
          if(sigma.cand > 0){
            if(length(lam0)==t){#if lam0 is year specific
              lamd.cand[,,l]<- lam0[l]*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
            }else{#fixed lam0
              lamd.cand[,,l]<- lam0*exp(-D[,,l]^2/(2*sigma.cand*sigma.cand))
            }
            pd.cand[,,l]=1-exp(-lamd.cand[,,l])
            ll.y.cand[,,l]= dbinom(y[,,l],K[l],pd.cand[,,l]*z[,l],log=TRUE) #only need to update this year
            ll.y.cand.t.sum[l]=sum(ll.y.cand[,,l])#just 1 year
            if(runif(1) < exp(ll.y.cand.t.sum[l] -ll.y.t.sum[l])){
              sigma[l]<- sigma.cand
              lamd[,,l]=lamd.cand[,,l]
              pd[,,l]=pd.cand[,,l]
              ll.y[,,l]=ll.y.cand[,,l]
              ll.y.t.sum[l]=ll.y.cand.t.sum[l]
            }
            ll.y.sum=sum(ll.y.t.sum)
          }
        }
      }else{#fixed sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D^2/(2*sigma.cand*sigma.cand))
          pd.cand=1-exp(-lamd.cand)
          for(l in 1:t){
            ll.y.cand[,,l]= dbinom(y[,,l],K[l],pd.cand[,,l]*z[,l],log=TRUE)
          }
          ll.y.cand.sum=sum(ll.y.cand)
          if(runif(1) < exp(ll.y.cand.sum - ll.y.sum)){
            sigma<- sigma.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
            ll.y.sum=ll.y.cand.sum
          }
        }
      }
      ##Update z[,1]
      if(t>3){
        upz=which(!(z[,1]==0&z[,2]==0&rowSums(z[,3:t]>0))&known.matrix[,1]==0)
      }else if(t==3){
        upz=which(!(z[,1]==0&z[,2]==0&z[,3]>0)&known.matrix[,1]==0)
      }else{
        upz=known.matrix[,1]==0
      }
      for(i in upz){
        z1.curr <- z[,1]
        z1.cand <- 1-z1.curr
        gamma.prime.cand <- gamma.prime
        z1.tmp <- z1.curr
        z1.tmp[i] <- z1.cand[i]
        a1.tmp <- a[,1]
        a1.tmp[i] <- 1-z1.cand[i]
        gamma.prime.cand[1] <- sum(z1.tmp)*gamma[1] / sum(a1.tmp)
        if(gamma.prime.cand[1] > 1) { # E(Recruits) must be < nAvailable
          warning("Rejected z due to low M")
          next
        }
        ll.z.cand[i,1] <- dbinom(z1.cand[i], 1, psi, log=TRUE)
        Ez.cand[i,1] <- z1.cand[i]*phi[1] + (1-z1.cand[i])*gamma.prime.cand[1]
        ll.z.cand[i,2] <- dbinom(z[i,2], 1, Ez.cand[i,1], log=TRUE)
        ll.y.cand[i,,1]=dbinom(y[i,,1],K[1],pd[i,,1]*z1.cand[i],log=TRUE)
        if(runif(1) < exp((sum(ll.y.cand[i,,1])+ sum(ll.z.cand[i,1:2]))-(sum(ll.y[i,,1])+sum(ll.z[i,1:2])) )) {#z1 and z2 matter
          ll.y[i,,1] <- ll.y.cand[i,,1]
          ll.z[i,1:2] <- ll.z.cand[i,1:2]
          Ez[i,1] <- Ez.cand[i,1]
          z[i,1] <- z1.cand[i]
          gamma.prime[1] <- gamma.prime.cand[1]
          if(z1.cand[i]==1){
            a[i,]=0#if caught occ 1, never available to recruit
          }else if (sum(z[i,])==0){#if never caught, turn availability all on
            a[i,]=1
          }else {#if not caught on occ1, but caught later, available to recruit in occ2
            a[i,1]=1
          }
        }
      }
      N[1]=sum(z[,1])

      #Update z[.2:t], gamma and phi
      # NOTE: Can't propose a 1 if next one is zero and subsequent is 1
      # Can't propose a 0 if 1 in known.matrix
      if(length(gamma)==1){
        gammause=c(gamma,rep(gamma,t-1))
      }else{
        gammause=gamma
      }
      if(length(phi)==1){
        phiuse=rep(phi,t-1)
      }else{
        phiuse=phi
      }
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
            rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)
            rem2=which(z[,l]==0&z[,l+1]==0&z[,t]>0)
          }else if(l==(t-1)){
            rem=which(rowSums(z[,1:(l-1)])>0&z[,t]>0)
            rem2=integer()
          }else if (l==t){
            rem=integer()
            rem2=integer()
          }else{
            rem=which(rowSums(z[,1:(l-1)])>0&rowSums(z[,(l+1):t])>0)
            rem2=which(z[,l]==0&z[,l+1]==0&rowSums(z[,(l+2):t])>0)
          }
        }else if(t==3){
          if(l==2){
            rem=which(z[,1]>0&z[,t]>0)
            rem2=integer()
          }else{#l==3
            rem=integer()
            rem2=integer()
          }
        }else{#t==2
          rem=integer()
          rem2=integer()
        }
        rem3=which( (z[,l-1]==0) & (a[,l-1]==0)) #dead guys
        remall=c(rem,rem2,rem3)
        upz=upz[!upz%in%remall]#We can change these guys
        navail=length(upz)
        if(navail<1){
          next
        }
        if(navail < proppars$propz[l-1]) {
          propz=navail
          warning("M isn't big enough")
        }else{
          propz=proppars$propz[l-1]
        }
        swapz=upz[sample.int(navail, propz)]
        pr.zt <- z[,l]
        pr.zt[swapz] <- Ez[swapz,l-1]
        zt.cand <- rbinom(M, 1, pr.zt)
        if(all(zt.cand == z[,l]))
          next
        at.cand=1*(a[,l-1]==1&zt.cand==0) #who was available on last occasion and not proposed to be captured?
        prop.probs <- sum(dbinom(zt.cand[swapz], 1, pr.zt[swapz], log=TRUE))
        back.probs <- sum(dbinom(z[swapz,t], 1, pr.zt[swapz], log=TRUE))
        ll.y.cand[swapz,,l] <- dbinom(y[swapz,,l], K[l],pd[swapz,,l]*zt.cand[swapz],log=TRUE)
        ll.z.cand[,l] <- dbinom(zt.cand, 1, Ez[,l-1], log=TRUE) ## Don't subset z
        prior.z <- sum(ll.z[,l])
        prior.z.cand <- sum(ll.z.cand[,l])
        if(l<t){ ## NOTE: Don't subset with swapz
          gamma.prime.cand[l] <- sum(zt.cand)*gammause[l-1] / sum(at.cand)
          if(gamma.prime.cand[l] > 1)
            next
          Ez.cand[,l] <- zt.cand*phiuse[l-1] + at.cand*gamma.prime.cand[l]
          ll.z.cand[,l+1] <- dbinom(z[,l+1], 1, Ez.cand[,l], log=TRUE)
          prior.z <- prior.z + sum(ll.z[,l+1])
          prior.z.cand <- prior.z.cand + sum(ll.z.cand[,l+1])
        }
        if(runif(1) < exp((sum(ll.y.cand[swapz,,l]) + prior.z.cand + back.probs) - (sum(ll.y[swapz,,l]) +prior.z + prop.probs) )) {
          ll.y[swapz,,l] <- ll.y.cand[swapz,,l]
          ll.z[,l] <- ll.z.cand[,l] ## NOTE: Don't subset with swapz
          z[,l] <- zt.cand
          a[,l] <- at.cand
          if(l < t) {
            ll.z[,l+1] <- ll.z.cand[,l+1]
            Ez[,l] <- Ez.cand[,l]
            if(t>3){#more a houskeeping if t>3
              #turn off availability if you died.
              if(l==2){
                dead=rowSums(z[swapz,l:t])==0&z[swapz,1:(l-1)]>0 #all 0 l and later but not all zero before
              }else if(l==t){
                dead=z[swapz,t]==0&rowSums(z[swapz,1:(l-1)])>0
              }else{
                dead=rowSums(z[swapz,l:t])==0&rowSums(z[swapz,1:(l-1)])>0
              }
              a[swapz[dead],l:t]=0
              #turn off availability for (l+1):t if you enter population
              a[swapz[z[swapz,l]==1],(l+1):t]=0
            }
            gamma.prime[l] <- gamma.prime.cand[l]
          }
        }
        N[l] <- sum(z[,l])
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
          gamma.prime.cand[l-1] <- N[l-1]*gamma.cand / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=!gamma.cand.ok
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
        }
        if(gamma.cand>0 & gamma.cand.ok) {
          #Only update ll.z for a=1 cases
          for(l in 2:t) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
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
          gamma.cand <- rnorm(1, gamma[l-1], proppars$gamma)
          gamma.cand.ok <- TRUE
          gamma.prime.cand[l-1] <- N[l-1]*gamma.cand / sum(a[,l-1])
          if(gamma.prime.cand[l-1] > 1){ ## Note don't break loop b/c ll.z needs updating because phi changed
            gamma.cand.ok=!gamma.cand.ok
          }
          Ez[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime[l-1]
          ll.z[,l] <- dbinom(z[,l], 1, Ez[,l-1], log=TRUE)
          if(gamma.cand>0 & gamma.cand.ok) {
            Ez.cand[,l-1] <- z[,l-1]*phiuse[l-1] + a[,l-1]*gamma.prime.cand[l-1]
            ll.z.cand[a[,l-1]==1,l] <- dbinom(z[a[,l-1]==1,l], 1, Ez.cand[a[,l-1]==1,l-1], log=TRUE)
            if(runif(1) < exp(sum(ll.z.cand[,l]) - sum(ll.z[,l]))) {
              gamma[l-1] <- gamma.cand
              gamma.prime[l-1] <- gamma.prime.cand[l-1]
              Ez[,l-1] <- Ez.cand[,l-1]
              ll.z[,l] <- ll.z.cand[,l]
            }
          }
        }
      }
      ## Now we have to update the activity centers
      if(metamu){
        #Update within year ACs
        for (j in 1:M) {
          for(l in 1:t){
            if(z[j,l]==0)
              next
            Scand <- c(rnorm(1, s2[j,l, 1], proppars$s2x), rnorm(1, s2[j,l, 2], proppars$s2y))
            if(useverts==FALSE){
              inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
            }else{
              inbox=inout(Scand,vertices)
            }
            if(inbox) {
              dtmp=sqrt((Scand[1] - X[[l]][, 1])^2 + (Scand[2] - X[[l]][, 2])^2)
              lamd.cand[j,1:nrow(X[[l]]),l]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
              pd.cand[j,,l]=1-exp(-lamd.cand[j,,l])
              ll.y.cand[j,,l] <- dbinom(y[j,,l], K[l], pd.cand[j,,l]*z[j,l], log=TRUE)
              ll.s2.cand<- dnorm(Scand[1],s1[j,1],sigma_t,log=TRUE)+dnorm(Scand[2],s1[j,2],sigma_t,log=TRUE)
              if(runif(1) < exp((sum(ll.y.cand[j,,l])+ll.s2.cand) -(sum(ll.y[j,,l])+ll.s2[j,l]))){
                s2[j,l,] <- Scand
                D[j,,l] <- dtmp
                lamd[j,,l] <- lamd.cand[j,,l]
                pd[j,,l]=pd.cand[j,,l]
                ll.y[j,,l]=ll.y.cand[j,,l]
                ll.s2[j,l]=ll.s2.cand
              }
            }
          }
        }
        #         #Update meta mus
        for (j in 1:M){
          if(sum(z[j,])==0)
            next
          Scand <- c(rnorm(1,s1[j,1],proppars$s1x), rnorm(1,s1[j,2],proppars$s1y))
          if(useverts==FALSE){
            inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
          }else{
            inbox=inout(Scand,vertices)
          }
          if(inbox){
            ll.s1=sum(ll.s2[j,])
            ll.s2cand<- dnorm(s2[j,,1],Scand[1],sigma_t,log=TRUE)+dnorm(s2[j,,2],Scand[2],sigma_t,log=TRUE)
            if (runif(1) < exp(sum(ll.s2cand) - sum(ll.s2[j,]))) {
              s1[j, ]=Scand
              ll.s2[j,]=ll.s2cand
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
        for (j in 1:M) {
          Scand <- c(rnorm(1, s1[j, 1], proppars$s2x), rnorm(1, s1[j, 2], proppars$s2y))
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
            lamd.cand[j,,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
            pd.cand[j,,]=1-exp(-lamd.cand[j,,])
            ll.y.cand[j,,]=ll.y[j,,]
            for(l in 1:t) {
              if(z[j,l]==0)
                next
              ll.y.cand[j,,l] <- dbinom(y[j,,l], K[l], pd.cand[j,,l]*z[j,l], log=TRUE)
            }
            if(runif(1) < exp(sum(ll.y.cand[j,,]) -sum(ll.y[j,,]))){
              s1[j, ] <- Scand
              D[j,, ] <- dtmp
              lamd[j,, ] <- lamd.cand[j,,]
              pd[j,,]=pd.cand[j,,]
              ll.y[j,,]=ll.y.cand[j,,]
            }
          }
        }
        for(l in 1:t){
          s2[,l,]=s1
        }
      }

      #Do we record output on this iteration?
      if(iter>(nburn)&iter%%nthin==0){
        s1xout[idx,]<- s1[,1]
        s1yout[idx,]<- s1[,2]
        zout[idx,,]<- z
        if(metamu){
          out[idx,]<- c(lam0,sigma ,gamma,phi,colSums(z),sigma_t)
          s2xout[idx,,]<- s2[,,1]
          s2yout[idx,,]<- s2[,,2]
        }else{
          out[idx,]<- c(lam0,sigma ,gamma,phi,N)
        }
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      if(metamu){
        list(out=out, s1xout=s1xout, s1yout=s1yout,s2xout=s2xout, s2yout=s2yout, zout=zout)
      }else{
        list(out=out, s1xout=s1xout, s1yout=s1yout, zout=zout)
      }
    }else{
      list(out=out)
    }
  }

