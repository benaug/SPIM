SCRmcmc <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=NA,obstype="poisson",
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           storeLatent=TRUE){
    ###
    library(abind)
    y<-data$y
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- data$K
    buff<- data$buff
    xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
    ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
    n=dim(y)[1]
    #Reduce to 2D data
    if(length(dim(y))==3){
      y=apply(y,c(1,2),sum)
    }
    
    #If using polygon state space
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
      xlim=c(min(vertices[,1]),max(vertices[,1]))
      ylim=c(min(vertices[,2]),max(vertices[,2]))
    }else if("buff"%in%names(data)){
      buff<- data$buff
      xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
      ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
      vertices=cbind(xlim,ylim)
      useverts=FALSE
    }else{
      stop("user must supply either 'buff' or 'vertices' in data object")
    }
    #Trap operation
    if("tf"%in%names(data)){
      tf=data$tf
      if(length(tf)!=nrow(X)){
        stop("If using a trap operation file, must input vector of length nrow(X).")
      }
 
    }else{
      tf=rep(K,J)
    }
    tf=matrix(rep(tf,M),ncol=J,nrow=M,byrow=TRUE)

    ##pull out initial values
    psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    known.vector=c(rep(1,n),rep(0,M-n))
    y=abind(y,matrix(0,nrow=M-n,ncol=J),along=1)
    #Initialize z
    z=1*(apply(y,1,sum)>0)
    add=M*(0.5-sum(z)/M)
    if(add>0){
      z[sample(which(z==0),add)]=1 #switch some uncaptured z's to 1.
    }
    
    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y[i,]>0,1:2],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s[i,]<- trps
      }
    }
    
    D=e2dist(s, X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y=array(0,dim=c(M,J))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      ll.y=dbinom(y,tf,pd*z,log=TRUE)
    }else if(obstype=="poisson"){
      ll.y=dpois(y,K*lamd*z,log=TRUE)
    }
    ll.y.cand=ll.y
    
    # some objects to hold the MCMC output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","psi"))
    if(storeLatent){
      sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
    }
    idx=1 #for storing output not recorded every iteration
    
    for(iter in 1:niter){
      #Update lam0
      if(obstype=="bernoulli"){
        llysum=sum(ll.y)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
          pd.cand=1-exp(-lamd.cand)
          ll.y.cand= dbinom(y,tf,pd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
            llysum=llycandsum
          }
        }
        #Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
          pd.cand=1-exp(-lamd.cand)
          ll.y.cand= dbinom(y,tf,pd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            sigma<- sigma.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
          }
        }
      }else{#poisson
        #Update lam0
        llysum=sum(ll.y)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
          ll.y.cand= dpois(y,tf*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
            llysum=llycandsum
          }
        }
        # Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
          ll.y.cand= dpois(y,tf*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            sigma<- sigma.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
          }
        }
      }
      
      ## probability of not being captured in a trap AT ALL
      if(obstype=="poisson"){
        pd=1-exp(-lamd)
      }
      pbar=(1-pd)^tf
      prob0<- exp(rowSums(log(pbar)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      if(obstype=="bernoulli"){
        ll.y= dbinom(y,tf,pd*z,log=TRUE)
      }else{
        ll.y= dpois(y,tf*lamd*z,log=TRUE)
      }
      psi=rbeta(1,1+sum(z),1+M-sum(z))
      
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
          lamd.cand[i,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          if(obstype=="bernoulli"){
            pd.cand[i,]=1-exp(-lamd.cand[i,])
            ll.y.cand[i,]= dbinom(y[i,],tf[i,],pd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              pd[i,]=pd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }else{#poisson
            ll.y.cand[i,]= dpois(y[i,],tf[i,]*lamd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(storeLatent){
          sxout[idx,]<- s[,1]
          syout[idx,]<- s[,2]
          zout[idx,]<- z
        }
        out[idx,]<- c(lam0,sigma ,sum(z),psi)
        idx=idx+1
      }
    }  # end of MCMC algorithm
    
    if(storeLatent){
      list(out=out, sxout=sxout, syout=syout, zout=zout)
    }else{ 
      list(out=out)
    }
  }

