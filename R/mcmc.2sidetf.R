mcmc.2sidetf <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,swap=10,swap.tol=1,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    ###
    library(abind)
    both<-data$both
    left<-data$left
    right<-data$right
    right<-data$right
    if(dim(right)[1]>dim(left)[1]){
      storeL=left[,2,,]
      storeR=right[,3,,]
      dimL=dim(left)
      left=array(0,dim=dim(right))
      right=array(0,dim=dimL)
      left[,2,,]=storeR
      right[,3,,]=storeL
      warning("Right side data set larger than left so I switched them for convenience")
    }
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(left)[3]
    IDknown<- data$IDknown
    Nfixed=length(IDknown)
    nleft<-dim(left)[1]-Nfixed
    nright<-dim(right)[1]-Nfixed
    nright<-dim(right)[1]-Nfixed
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
      ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    #trap history
    tf=t(data$tf)
    tfnew=abind(tf,tf,along=0)
    for(i in 3:M){
      tfnew=abind(tfnew,tf,along=1)
    }
    tf=tfnew
    #Figure out what needs to be updated
    lam01<- inits$lam01
    lam02<- inits$lam02
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

    ##pull out initial values
    psi<- inits$psi
    sigma<- inits$sigma
    if(uplam02){
      lam02<- inits$lam02
    }else{
      lam02=0
      if(!is.na(inits$lam02)){
        warning("User specified init for lam02 with no double traps. Setting to NA.")
      }
    }

    #augment both data
    both<- abind(both,array(0, dim=c( M-dim(both)[1],3,K, J)), along=1)
    left<- abind(left,array(0, dim=c( M-dim(left)[1],3, K, J)), along=1)
    right<-abind(right,array(0,dim=c( M-dim(right)[1],3,K,J)), along=1)

    #sort to minimize distance between initial matches. Skip if no single sides. Need to fix if only lefts or only rights >Nfixed
    if(nleft>0|nright>0){
      IDs<- LRmatch(M=M,left=left, nleft=nleft, right=right, nright=nright, X, Nfixed=Nfixed)
      #Add unused augmented indivuals back in
      notusedL<- (1:M)[is.na(match(1:M,IDs$ID_L))]
      ID_L<-c(IDs$ID_L,notusedL)
      notusedR<- (1:M)[is.na(match(1:M,IDs$ID_R))]
      ID_R<-c(IDs$ID_R,notusedR)
    }else{
      ID_R=ID_L=1:M
    }

    #Make initial complete data set
    tmpdata<- both + left[order(ID_L),,,] + right[order(ID_R),,,]
    tmpdata<- apply(tmpdata,c(1,4),sum)
    z=1*(apply(tmpdata,1,sum)>0)
    z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?

    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
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
    idx=which(rowSums(tmpdata)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[tmpdata[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }

    known.vector<- c( rep(1,Nfixed), rep(0, M-Nfixed) )

    #Bernoulli Likelihood function
    func<- function(lamd1,lamd2,y.both,y.left,y.right,tf,z){
      #convert lamd to pd (gaussian hazard model)
      pd1=1-exp(-lamd1)
      pd2=1-exp(-lamd2)
      ones=tf==1
      twos=tf==2
      part_both <- dbinom(y.both*twos,tf==2,pd2*twos,log=TRUE)
      pd1b=2*pd1-pd1*pd1
      part_left <-  dbinom(y.left*ones,ones,pd1*ones,log=TRUE)+dbinom(y.left*twos,twos,pd1b*twos,log=TRUE)
      part_right <-  dbinom(y.right*ones,ones,pd1*ones,log=TRUE)+dbinom(y.right*twos,twos,pd1b*twos,log=TRUE)
      v<- part_both+part_left+part_right
      #If data is M x K or 2 x K
      if(length(dim(y.both))==3){
        v[z==0,,]<- 0
      }else{
        v<- v*z
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
    sxout<- syout<- zout<- ID_Lout<- ID_Rout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration
    #3D data when using trap file. Need to know when double traps go to single traps
    y.both<- apply(both[,1,,], c(1,2,3), sum)
    y.left<- apply(left[order(ID_L),2,,], c(1,2,3), sum)
    y.right<- apply(right[order(ID_R),3,,], c(1,2,3), sum)
    D<- e2dist(s, X)
    Dnew=abind(D,D,along=0)
    for(k in 3:K){
      Dnew=abind(Dnew,D,along=1)
    }
    D=aperm(Dnew,c(2,1,3))
    lamd1<- lam01*exp(-D*D/(2*sigma*sigma))
    lamd2<- lam02*exp(-D*D/(2*sigma*sigma))
    zero.guys<- apply(y.both+y.left + y.right ,1,sum) == 0
    lik.curr<-  sum( func(lamd1,lamd2,y.both, y.left,y.right,tf,z) )

    for(i in 1:niter){
      #Update lam01
      if(uplam01){
        lik.curr<-  sum( func(lamd1,lamd2,y.both, y.left,y.right,tf,z) )
        lam01.cand<- rnorm(1,lam01,proppars$lam01)
        if(lam01.cand > 0){
          lamd1.cand<- lam01.cand*exp(-D*D/(2*sigma*sigma))
          lik.new<-  sum( func(lamd1.cand,lamd2,y.both,y.left,y.right,tf,z) )
          if(runif(1) < exp(lik.new -lik.curr)){
            lam01<- lam01.cand
            lamd1=lamd1.cand
            lik.curr<- lik.new
          }
        }
      }
      #Update lam02
      if(uplam02){
        lam02.cand<- rnorm(1,lam02,proppars$lam02)
        if(lam02.cand > 0){
          lamd2.cand<- lam02.cand*exp(-D*D/(2*sigma*sigma))
          lik.new<-  sum( func(lamd1,lamd2.cand,y.both,y.left,y.right,tf,z) )
          if(runif(1) < exp(lik.new -lik.curr)){
            lam02<- lam02.cand
            lamd2=lamd2.cand
            lik.curr<- lik.new
          }
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamd1.cand<- lam01*exp(-D*D/(2*sigma.cand*sigma.cand))
        lamd2.cand<- lam02*exp(-D*D/(2*sigma.cand*sigma.cand))
        lik.new<-  sum( func(lamd1.cand,lamd2.cand,y.both,y.left,y.right,tf,z) )
        if(runif(1) < exp(lik.new -lik.curr)){
          sigma<- sigma.cand
          lamd1=lamd1.cand
          lamd2=lamd2.cand
          lik.curr<- lik.new
        }
      }
      #
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
        candmap[z==0,1]=NA #Don't swap in z=0 guys
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
          newID<-ID_L
          newID[guy1]<- s.swap.in
          newID[guy2]<- s.swap.out

          ## recompute 'data' and compute likelihood components for swapped guys only
          y.left.tmp<- apply( left[order(newID),2,,], c(1,2,3), sum)#Only y.left changed
          swapped=c(s.swap.out,s.swap.in)
          llswappre<- sum(func(lamd1[swapped,,],lamd2[swapped,,],y.both[swapped,,],y.left[swapped,,],y.right[swapped,,],tf[swapped,,],z[swapped]))
          llswappost<- sum( func(lamd1[swapped,,],lamd2[swapped,,],y.both[swapped,,],y.left.tmp[swapped,,],y.right[swapped,,],tf[swapped,,],z[swapped]))
          lldiff=llswappost-llswappre
          #MH step
          if(runif(1)<exp(lldiff)*(jump.back/jump.probability) ){
            lik.curr=lik.curr+lldiff #update likelihood
            y.left <- y.left.tmp #update left data
            ID_L<-newID #update data order
            map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
            zero.guys<- apply(y.both+y.left + y.right ,1,sum) == 0
            if( sum(y.both[guy1,,] + y.left.tmp[guy1,,] )>0 & z[guy1]==0){
              cat("error on guy1",fill=TRUE)
              browser()
            }
            if( sum(y.both[guy2,,] + y.left.tmp[guy2,,] )>0 & z[guy2]==0){
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
          newID<-ID_R
          newID[guy1]<- s.swap.in
          newID[guy2]<- s.swap.out

          ## recompute 'data' and compute likelihood components for swapped guys only
          y.right.tmp<- apply( right[order(newID),3,,], c(1,2,3), sum) #Only y.right changed
          swapped=c(s.swap.out,s.swap.in)
          llswappre<- sum(func(lamd1[swapped,,],lamd2[swapped,,],y.both[swapped,,],y.left[swapped,,],y.right[swapped,,],tf[swapped,,],z[swapped]))
          llswappost<- sum( func(lamd1[swapped,,],lamd2[swapped,,],y.both[swapped,,],y.left[swapped,,],y.right.tmp[swapped,,],tf[swapped,,],z[swapped]))
          lldiff=llswappost-llswappre
          #MH step
          if(runif(1)<exp(lldiff)*(jump.back/jump.probability) ){
            lik.curr=lik.curr+lldiff #update likelihood
            y.right <- y.right.tmp #update right data
            ID_R<-newID #update data order
            map[c(guy1,guy2),2]=c(s.swap.in,s.swap.out)
            zero.guys<- apply(y.both+y.left + y.right ,1,sum) == 0
            if( sum(y.both[guy1,,] + y.right.tmp[guy1,,] )>0 & z[guy1]==0){
              cat("error on guy1",fill=TRUE)
              browser()
            }
            if( sum(y.both[guy2,,] + y.right.tmp[guy2,,] )>0 & z[guy2]==0){
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
      #pbar0<- (1 - lam01*exp(-D*D/(2*sigma*sigma)))^2*K * (1 - lam02*exp(-D*D/(2*sigma*sigma)))^K
      trapno=matrix(rep(X[,3],M),nrow=M,byrow=TRUE) #trap number multiplier for left and right captures
      ones=trapno==1
      twos=trapno==2
      pd1=1-exp(-lamd1)
      pd2=1-exp(-lamd2)
      pd1b=(2*pd1-pd1*pd1)
      pd1traps=pd1*(tf==1)+pd1b*(tf==2)
      pd2traps=pd2*(tf==2)
      pbar1=apply(1-pd1traps,1,prod)^2
      pbar2=apply(1-pd2traps,1,prod)
      prob0=pbar1*pbar2
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[zero.guys & known.vector==0]<- rbinom(sum(zero.guys & (known.vector ==0) ), 1, fc[zero.guys & (known.vector==0) ])
      lik.curr<-  sum( func(lamd1,lamd2,y.both, y.left,y.right,tf,z) )
      psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

      ## Now we have to update the activity centers
      for (j in 1:M) {
        Scand <- c(rnorm(1, s[j, 1], proppars$sx), rnorm(1, s[j, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          dtmp=matrix(rep(dtmp,K),ncol=J,nrow=K,byrow=T)
          lamd1.thisj<- lam01*exp(-D[j,,]*D[j,,]/(2*sigma*sigma))
          lamd2.thisj<- lam02*exp(-D[j,,]*D[j,,]/(2*sigma*sigma))
          lamd1.cand<- lam01*exp(-dtmp*dtmp/(2*sigma*sigma))
          lamd2.cand<- lam02*exp(-dtmp*dtmp/(2*sigma*sigma))
          llS<- sum(func(lamd1.thisj,lamd2.thisj,y.both[j,,],y.left[j,,],y.right[j,,],tf[j,,],z[j]))
          llcand<- sum(func(lamd1.cand,lamd2.cand,y.both[j,,],y.left[j,,],y.right[j,,],tf[j,,],z[j]))

          if (runif(1) < exp(llcand - llS)) {
            s[j, ] <- Scand
            D[j,, ] <- dtmp
            lamd1[j,, ] <- lamd1.cand
            lamd2[j,,] <- lamd2.cand
          }
        }
      }
      #Do we record output on this iteration?
      if(i>nburn&i%%nthin==0){
        sxout[idx,]<- s[,1]
        syout[idx,]<- s[,2]
        zout[idx,]<- z
        ID_Lout[idx,]<-ID_L
        ID_Rout[idx,]<-ID_R
        out[idx,]<- c(lam01,lam02,sigma ,sum(z))
        idx=idx+1
      }
    }  # end of MCMC algorithm
    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout, ID_Lout=ID_Lout,ID_Rout=ID_Rout)
    }else{
      list(out=out)
    }
  }
