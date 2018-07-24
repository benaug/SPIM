mcmc.2sideR <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,
           swap=10,swap.tol=1,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),storeLatent=TRUE){
    library(abind)
    y.both<-data$both
    y.left.obs<-data$left
    y.right.obs<-data$right
    if(length(dim(y.both))==3){
      y.both=apply(y.both,c(1,2),sum)
      y.left.obs=apply(y.left.obs,c(1,2),sum)
      y.right.obs=apply(y.right.obs,c(1,2),sum)
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
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- data$K
    if("tf"%in%names(data)){
      tf=data$tf
      if(is.matrix(tf)){
        stop("This is the wrong function for a 2D trap file. Tell Ben something went wrong")
      }
    }else{
      tf=rep(K,J)
    }
    tf=matrix(rep(tf,M),nrow=M,byrow=TRUE)
    IDknown<- data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(y.left.obs)[1]-Nfixed
    nright<-dim(y.right.obs)[1]-Nfixed
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
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
    uplam01=uplam02=upIDs=TRUE
    if(lam01==0){
      uplam01=FALSE
      upIDs=FALSE
    }
    if(all(X[,3]==1)|lam02==0){
      uplam02=FALSE
    }
    if(upIDs==TRUE&(nleft==0&nright==0)){  #not sure what will happen if we only have left only or right only
      upIDs=FALSE
    }

    #augment 3 data sets
    y.both<- abind(y.both,array(0, dim=c( M-dim(y.both)[1],J)), along=1)
    y.left.obs<- abind(y.left.obs,array(0, dim=c( M-dim(y.left.obs)[1],J)), along=1)
    y.right.obs<-abind(y.right.obs,array(0,dim=c( M-dim(y.right.obs)[1],J)), along=1)

    #sort to minimize distance between initial matches. Skip if no single sides.
    if(nleft>0|nright>0){
      IDs<- LRmatch(M=M,left=y.left.obs, nleft=nleft, right=y.right.obs, nright=nright, X, Nfixed=Nfixed)
      #Add unused augmented indivuals back in
      notusedL<- (1:M)[is.na(match(1:M,IDs$ID_L))]
      ID_L<-c(IDs$ID_L,notusedL)
      notusedR<- (1:M)[is.na(match(1:M,IDs$ID_R))]
      ID_R<-c(IDs$ID_R,notusedR)
    }else{
      ID_R=ID_L=1:M
    }
    
    #reoder left and right possibly true data sets
    y.left.true<- y.left.obs[order(ID_L),]
    y.right.true<- y.right.obs[order(ID_R),]
    
    #Make initial complete data set
    tmpdata<- y.both + y.left.true + y.right.true
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
        inside[i]=Rcpp::inout(s[i,],vertices)
      }
      idx=which(inside==FALSE)
      if(length(idx)>0){
        for(i in 1:length(idx)){
          while(inside[idx[i]]==FALSE){
            s[idx[i],]=c(runif(1,xlim[1],xlim[2]), runif(1,ylim[1],ylim[2]))
            inside[idx[i]]=Rcpp::inout(s[idx[i],],vertices)
          }
        }
      }
    }
    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam01","lam02","sigma","N"))
    if(storeLatent){
      sxout<- syout<- zout<- ID_Lout<- ID_Rout<-matrix(NA,nrow=nstore,ncol=M)
    }
    idx=1 #for storing output not recorded every iteration
    zero.guys<- apply(y.both+y.left.true + y.right.true ,1,sum) == 0
    trapno=matrix(rep(X[,3],2),nrow=M,ncol=J,byrow=TRUE) #trap number multiplier for left and right captures.
    ones=trapno==1
    twos=trapno==2
    
    D<- e2dist(s, X)
    lamd1<- lam01*exp(-D*D/(2*sigma*sigma))
    lamd2<- lam02*exp(-D*D/(2*sigma*sigma))
    pd1=1-exp(-lamd1)
    pd1b=ones*pd1+twos*(2*pd1-pd1*pd1)
    pd2=1-exp(-lamd2)
    lamd1.cand=lamd1
    lamd2.cand=lamd2
    pd1.cand=pd1
    pd2.cand=pd2
    pd1b.cand=pd1b

    ll.y.both <- dbinom(y.both,tf,z*pd2*twos,log=TRUE)
    ll.y.left <-  dbinom(y.left.true,tf,z*pd1b,log=TRUE)
    ll.y.right <-  dbinom(y.right.true,tf,z*pd1b,log=TRUE)
    ll.y.both.cand=ll.y.both
    ll.y.left.cand=ll.y.left
    ll.y.right.cand=ll.y.right
    if(!is.finite(sum(ll.y.both)))stop("Both side likelihood not finite. Make sure all camera stations recording both side captures have 2 cameras. Then try changing lam02 or sigma inits.")
    if(!is.finite(sum(ll.y.left)))stop("Left side likelihood not finite. Try changing lam01 or sigma inits.")
    if(!is.finite(sum(ll.y.right)))stop("right side likelihood not finite. Try changing lam01 or sigma inits.")
    
    for(iter in 1:niter){
      #Update lam01
      lly1sum=sum(ll.y.left)+sum(ll.y.right)
      if(uplam01){
        lam01.cand<- rnorm(1,lam01,proppars$lam01)
        if(lam01.cand > 0){
          lamd1.cand<- lam01.cand*exp(-D*D/(2*sigma*sigma))
          pd1.cand=1-exp(-lamd1.cand)
          pd1b.cand=ones*pd1.cand+twos*(2*pd1.cand-pd1.cand*pd1.cand)
          ll.y.left.cand=dbinom(y.left.true,tf,z*pd1b.cand,log=TRUE)
          ll.y.right.cand=dbinom(y.right.true,tf,z*pd1b.cand,log=TRUE)
          lly1candsum=sum(ll.y.left.cand)+sum(ll.y.right.cand)
          if(runif(1) < exp(lly1candsum - lly1sum)){
            lam01<- lam01.cand
            lamd1=lamd1.cand
            pd1=pd1.cand
            pd1b=pd1b.cand
            ll.y.left=ll.y.left.cand
            ll.y.right=ll.y.right.cand
            lly1sum=lly1candsum
          }
        }
      }
      #Update lam02
      lly2sum=sum(ll.y.both)
      if(uplam02){
        lam02.cand<- rnorm(1,lam02,proppars$lam02)
        if(lam02.cand > 0){
          lamd2.cand<- lam02.cand*exp(-D*D/(2*sigma*sigma))
          pd2.cand=1-exp(-lamd2.cand)
          ll.y.both.cand <- dbinom(y.both,tf,z*pd2.cand*twos,log=TRUE)
          lly2candsum=sum(ll.y.both.cand)
          if(runif(1) < exp(lly2candsum-lly2sum)){
            lam02=lam02.cand
            lamd2=lamd2.cand
            pd2=pd2.cand
            ll.y.both=ll.y.both.cand
            lly2sum=lly2candsum
          }
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamd1.cand<- lam01*exp(-D*D/(2*sigma.cand*sigma.cand))
        lamd2.cand<- lam02*exp(-D*D/(2*sigma.cand*sigma.cand))
        pd1.cand=1-exp(-lamd1.cand)
        pd1b.cand=ones*pd1.cand+twos*(2*pd1.cand-pd1.cand*pd1.cand)
        pd2.cand=1-exp(-lamd2.cand)
        ll.y.left.cand=dbinom(y.left.true,tf,z*pd1b.cand,log=TRUE)
        ll.y.right.cand=dbinom(y.right.true,tf,z*pd1b.cand,log=TRUE)
        lly1candsum=sum(ll.y.left.cand)+sum(ll.y.right.cand)
        ll.y.both.cand <- dbinom(y.both,tf,z*pd2.cand*twos,log=TRUE)
        lly2candsum=sum(ll.y.both.cand)
        if(runif(1) < exp((lly1candsum+lly2candsum) -(lly1sum+lly2sum))){
          sigma<- sigma.cand
          lamd1=lamd1.cand
          lamd2=lamd2.cand
          pd1=pd1.cand
          pd1b=pd1b.cand
          pd2=pd2.cand
          ll.y.both=ll.y.both.cand
          ll.y.left=ll.y.left.cand
          ll.y.right=ll.y.right.cand
        }
      }

      ## Update latent ID variables
      ## Candidate: swap IDs of one guy with some other guy. The two guys to be swapped are
      ## chosen randomly from the z=1 guys
      if(any(z[!zero.guys]==0)){
        cat("some z = 0!",fill=TRUE)
        browser()
      }
      #Swap left guys first
      if(nleft>0){
        # User inputs how many swaps to make on each iteration my specifying "swap"
        #map lefts to boths
        #candmap used to remove disallowed candidates.NA any lefts that map back to z=0 indices and boths that are z=0 indices
        map=cbind(1:M,ID_L)
        candmap=map
        candmap[1:Nfixed,]=NA #Don't swap out IDknown guys
        candmap[z==0,1]=NA #Don't swap in z=0 guys.
        candmap[candmap[,2]%in%which(z==0),2]=NA #Don't choose a guy1 that will make you swap in a z=0 guy
        #These are the guys that can come out and go in
        OUTcands=which(!is.na(candmap[,2]))
        INcands=which(!is.na(candmap[,1]))
        for(up in 1:swap){
          ### In this code here "guy1" and "guy2" are indexing "left-side encounter histories" whereas
          ### "s.swap.in" and "s.swap.out" are indexing (both-side guys, s, z)
          guy1<- sample(OUTcands, 1)
          s.swap.out<- map[guy1,2]
          # to find candidates for swapping look in vicinity of it's current both side membership
          dv<-  sqrt( (s[s.swap.out,1]- s[INcands,1])^2 + (s[s.swap.out,2] - s[INcands,2])^2 )
          ## This is a list of both side activity centers that could be assinged to the right side guy
          ## under consideration based on how far away they are.  swap.tol is the spatial tolerance
          # if no one around, code swaps guy 1 with guy 1 (no swap, but does calculations)
          possible<- INcands[dv < swap.tol]
          jump.probability<- 1/length(possible)  # h(theta*|theta)
          # this is a particular value of s to swap for guy1.id
          if(length(possible)>1){
            s.swap.in <-  sample( possible, 1)
          }
          if(length(possible) ==1){
            s.swap.in <- possible
          }
          if(length(possible)==0) next #no guys close enough to swap

          #compute  h(theta|theta*)
          trash<-   sqrt( (s[s.swap.in,1]- s[INcands,1])^2 + (s[s.swap.in,2] - s[INcands,2])^2  )
          trash<-  INcands[trash < swap.tol]
          jump.back<-  1/length(trash)

          ##  Which left encounter history is currently associated with both guy s.swap.in?
          guy2<- which(map[,2]==s.swap.in)
          if(guy1==guy2) next #Same guy selected
          newID<-ID_L
          newID[guy1]<- s.swap.in
          newID[guy2]<- s.swap.out

          ## recompute 'data' and compute likelihood components for swapped guys only
          y.left.tmp<- y.left.obs[order(newID),]#Reorder observed left side data set
          swapped=c(s.swap.out,s.swap.in)
          ll.y.left.tmp=dbinom(y.left.tmp[swapped,],tf[swapped,],pd1b[swapped,],log=TRUE)
    
          #MH step
          if(runif(1)<exp(sum(ll.y.left.tmp)-sum(ll.y.left[swapped,]))*(jump.back/jump.probability) ){
            y.left.true=y.left.tmp #update left data
            ll.y.left[swapped,]=ll.y.left.tmp
            ID_L<-newID #update data order
            map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
            zero.guys<- apply(y.both+y.left.true + y.right.true ,1,sum) == 0
            if( sum(y.both[guy1,] + y.left.tmp[guy1,] )>0 & z[guy1]==0){
              cat("error on guy1",fill=TRUE)
              browser()
            }
            if( sum(y.both[guy2,] + y.left.tmp[guy2,] )>0 & z[guy2]==0){
              cat("error on guy2",fill=TRUE)
              browser()
            }
          }
        }
        #Do any captured guys have z==0?
        if(any(z[!zero.guys]==0)){
          cat("coming out of ID_L update: some z = 0!",fill=TRUE)
          browser()
        }
      }

      #Repeat for rights
      if(nright>0){
        #Reuse left side of maps
        map[,2]=candmap[,2]=ID_R
        candmap[1:Nfixed,2]=NA
        candmap[candmap[,2]%in%which(z==0),2]=NA
        OUTcands=which(!is.na(candmap[,2])) #INcands the same
        for(up in 1:swap){
          ### In this code here "guy1" and "guy2" are indexing "right-side encounter histories" whereas
          ### "s.swap.in" and "s.swap.out" are indexing (both-side guys, s, z)
          guy1<- sample(OUTcands, 1)
          s.swap.out<- map[guy1,2]
          # to find candidates for swapping look in vicinity of it's current both side membership
          dv<-  sqrt( (s[s.swap.out,1]- s[INcands,1])^2 + (s[s.swap.out,2] - s[INcands,2])^2 )
          ## This is a list of both side activity centers that could be assinged to the right side guy
          ## under consideration based on how far away they are.  swap.tol is the spatial tolerance
          possible<- INcands[dv < swap.tol]
          jump.probability<- 1/length(possible)  # h(theta*|theta)
          # this is a particular value of s to swap for guy1.id
          if(length(possible)>1){
            s.swap.in <-  sample( possible, 1)
          }
          if(length(possible) ==1){
            s.swap.in <- possible
          }
          if(length(possible)==0) next #no guys close enough to swap

          #compute  h(theta|theta*)
          trash<-   sqrt( (s[s.swap.in,1]- s[INcands,1])^2 + (s[s.swap.in,2] - s[INcands,2])^2  )
          trash<-  INcands[trash < swap.tol]
          jump.back<-  1/length(trash)

          ##  Which right encounter history is currently associated with both guy s.swap.in?
          guy2<- which(map[,2]==s.swap.in)
          if(guy1==guy2) next #Same guy selected
          newID<-ID_R
          newID[guy1]<- s.swap.in
          newID[guy2]<- s.swap.out
          
          ## recompute 'data' and compute likelihood components for swapped guys only
          y.right.tmp<- y.right.obs[order(newID),]#Reorder observed right side data set
          swapped=c(s.swap.out,s.swap.in)
          ll.y.right.tmp=dbinom(y.right.tmp[swapped,],tf[swapped,],pd1b[swapped,],log=TRUE)
          
          #MH step
          if(runif(1)<exp(sum(ll.y.right.tmp)-sum(ll.y.right[swapped,]))*(jump.back/jump.probability) ){
            y.right.true=y.right.tmp #update right data
            ll.y.right[swapped,]=ll.y.right.tmp
            ID_R<-newID #update data order
            map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
            zero.guys<- apply(y.both+y.left.true + y.right.true ,1,sum) == 0
            if( sum(y.both[guy1,] + y.right.tmp[guy1,] )>0 & z[guy1]==0){
              cat("error on guy1",fill=TRUE)
              browser()
            }
            if( sum(y.both[guy2,] + y.right.tmp[guy2,] )>0 & z[guy2]==0){
              cat("error on guy2",fill=TRUE)
              browser()
            }
          }
        }
        #Do any captured guys have z==0?
        if(any(z[!zero.guys]==0)){
          cat("coming out of ID_R update: some z = 0!",fill=TRUE)
          browser()
        }
      }
      #Update psi gibbs
      ## probability of not being captured in a trap AT ALL
      pd2b=pd2*twos
      pbar1=(1-pd1b)^(2*tf) #can be a left or a right
      pbar2=(1-pd2b)^tf
      pbar0=pbar1*pbar2
      prob0<- exp(rowSums(log(pbar0)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[zero.guys & known.vector==0]<- rbinom(sum(zero.guys & (known.vector ==0) ), 1, fc[zero.guys & (known.vector==0) ])
      #update observation model ll
      ll.y.both <- dbinom(y.both,tf,z*pd2b,log=TRUE)
      ll.y.left <-  dbinom(y.left.true,tf,z*pd1b,log=TRUE)
      ll.y.right <-  dbinom(y.right.true,tf,z*pd1b,log=TRUE)
      psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

      ## Now we have to update the activity centers
      for (i in 1:M) {
        Scand <- c(rnorm(1, s[i, 1], proppars$sx), rnorm(1, s[i, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamd1.cand[i,]=lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          lamd2.cand[i,]=lam02*exp(-dtmp*dtmp/(2*sigma*sigma))
          pd1.cand[i,]=1-exp(-lamd1.cand[i,])
          pd1b.cand[i,]=ones[i,]*pd1.cand[i,]+twos[i,]*(2*pd1.cand[i,]-pd1.cand[i,]*pd1.cand[i,])
          pd2.cand[i,]=1-exp(-lamd2.cand[i,])
          ll.y.both.cand[i,]=dbinom(y.both[i,],tf[i,],z[i]*pd2.cand[i,]*twos[i,],log=TRUE)
          ll.y.left.cand[i,]=dbinom(y.left.true[i,],tf[i,],z[i]*pd1b.cand[i,],log=TRUE)
          ll.y.right.cand[i,]=dbinom(y.right.true[i,],tf[i,],z[i]*pd1b.cand[i,],log=TRUE)
          llysum=(sum(ll.y.both[i,])+sum(ll.y.left[i,])+sum(ll.y.right[i,]))
          llycandsum=(sum(ll.y.both.cand[i,])+sum(ll.y.left.cand[i,])+sum(ll.y.right.cand[i,]))
          if (runif(1) < exp(llycandsum-llysum)) {
            s[i, ] <- Scand
            D[i, ] <- dtmp
            lamd1[i, ] <- lamd1.cand[i,]
            lamd2[i, ] <- lamd2.cand[i,]
            pd1[i,]=pd1.cand[i,]
            pd1b[i,]=pd1b.cand[i,]
            pd2[i,]=pd2.cand[i,]
            ll.y.both[i,]=ll.y.both.cand[i,]
            ll.y.left[i,]=ll.y.left.cand[i,]
            ll.y.right[i,]=ll.y.right.cand[i,]
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(storeLatent){
          sxout[idx,]<- s[,1]
          syout[idx,]<- s[,2]
          zout[idx,]<- z
          ID_Lout[idx,]<-ID_L
          ID_Rout[idx,]<-ID_R
        }
        out[idx,]<- c(lam01,lam02,sigma ,sum(z))
        idx=idx+1
      }
    }  # end of MCMC algorithm
    if(storeLatent==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout, ID_Lout=ID_Lout,ID_Rout=ID_Rout)
    }else{
      list(out=out)
    }
  }