SCRIPPmcmc <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2,beta0=0.05,beta1=0.05),cellArea,keepACs=TRUE){
    ###
    library(abind)
    y<-data$y
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(y)[2]
    n<- data$n
    grid=as.matrix(data$grid)
    npix=nrow(grid)
    ##pull out initial values
    beta0<- inits$beta0
    beta1<- inits$beta1
    lam0<- inits$lam0
    sigma<- inits$sigma

    #augment data
    y<- abind(y,array(0, dim=c( M-dim(y)[1],K, J)), along=1)
    known.vector=c(rep(1,data$n),rep(0,M-data$n))

    #Make initial complete data set
    y2D=apply(y,c(1,3),sum)

    #Optimize starting locations given where they are trapped.
    cell=sample(1:npix,M)#assign random locations
    s<- cbind(cell,grid[cell,1:2])
    idx=which(rowSums(y)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[y2D[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      loc=c(mean(trps[,1]),mean(trps[,2]))
      d=sqrt((loc[1]-grid[,1])^2+(loc[2]-grid[,2])^2)
      cell=which(d==min(d))
      if(length(cell)>1){
        cell=cell[1]
      }
      s[i,]=c(cell,grid[cell,1:2])
    }

    #Bernoulli Likelihood function
    func<- function(lamd,y,K,z,X){
      #convert lamd to pd (gaussian hazard model)
      pd=1-exp(-lamd)
      #If data is M x K
      if(is.matrix(y)){
        v <-  dbinom(y,K,pd,log=TRUE)
        v[z==0,]<- 0
      }else{
        #If data is 1 x K
        v <- dbinom(y,K,pd,log=TRUE)
        v<- v*z
      }
      v
    }

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=6)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","EN","beta0","beta1"))
    sxout<- syout<- scellout<-zout<-matrix(NA,nrow=nstore,ncol=M)
    idx=1 #for storing output not recorded every iteration

    D<- e2dist(s[,2:3], X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))

    mu=exp(beta0+beta1*1*(grid[,3]==2))*cellArea
    EN=sum(mu)
    psi <- EN/M
    if (psi > 1){
      stop("Bad initial values for beta0 or beta1. Or M is too low")
    }
    z=1*(apply(y2D,1,sum)>0)
    N=sum(z)
    N2=round(psi*EN)

    if(N<N2){
      z[sample(which(z==0),N2-N)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?
      N=sum(z)
    }

    for(i in 1:niter){
      #Update lam0
      lik.curr<-  sum( func(lamd,y2D,K,z,X) )
      lam0.cand<- rnorm(1,lam0,proppars$lam0)
      if(lam0.cand > 0){
        lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
        lik.new<-  sum( func(lamd.cand,y2D,K,z,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          lam0<- lam0.cand
          lamd=lamd.cand
          lik.curr<- lik.new
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
        lik.new<-  sum( func(lamd.cand,y2D,K,z,X) )
        if(runif(1) < exp(lik.new -lik.curr)){
          sigma<- sigma.cand
          lamd=lamd.cand
          lik.curr<- lik.new
        }
      }
      #update z
      ## probability of not being captured in a trap AT ALL
      pd=1-exp(-lamd)
      pbar=(1-pd)^K
      prob0<- exp(rowSums(log(pbar)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      lik.curr<-  sum( func(lamd,y,K,z,X) )
      N=sum(z)

      #Update beta0 and beta1
      beta0.cand <- rnorm(1, beta0, proppars$beta0)
      mu.cand=exp(beta0.cand+beta1*(grid[,3]==2))*cellArea
      EN.cand=sum(mu.cand)
      ll.beta <-   sum(((beta0 + beta1 * (grid[s[,1],3]==2)) - log(EN)) * z) + dbinom(N, M, EN/M, log = TRUE)
      #ll.beta <-   sum((beta0 + beta1 * (grid[s[,1],3]==2))* z)-EN + dbinom(N, M, EN/M, log = TRUE)
      if (EN.cand < M) {
        ll.beta.cand <- sum(((beta0.cand + beta1 * (grid[s[,1],3]==2)) - log(EN.cand)) * z) + dbinom(N, M, EN.cand/M,log = TRUE)
        #ll.beta.cand <- sum((beta0.cand + beta1 * (grid[s[,1],3]==2)) * z)- EN.cand + dbinom(N, M, EN.cand/M,log = TRUE)
        if (runif(1) < exp(ll.beta.cand - ll.beta)) {
          beta0 <- beta0.cand
          EN <- EN.cand
          ll.beta <- ll.beta.cand
        }
      }
      beta1.cand <- rnorm(1, beta1, proppars$beta1)
      mu.cand=exp(beta0+beta1.cand*(grid[,3]==2))*cellArea
      EN.cand=sum(mu.cand)
      if (EN.cand < M) {
        #ll.beta.cand <- sum((beta0 + beta1.cand * (grid[s[,1],3]==2)) * z) - EN.cand + dbinom(N, M, EN.cand/M,log = TRUE)
        ll.beta.cand <- sum(((beta0 + beta1.cand * (grid[s[,1],3]==2)) - log(EN.cand)) * z) + dbinom(N, M, EN.cand/M,log = TRUE)
        if (runif(1) < exp(ll.beta.cand - ll.beta)) {
          beta1 <- beta1.cand
          EN <- EN.cand
          ll.beta <- ll.beta.cand
        }
      }
      psi <- EN/M

      ## Now we have to update the activity centers
      for (j in 1:M) {
        Scand <- c(rnorm(1, s[j, 2], proppars$sx), rnorm(1, s[j, 3], proppars$sy))
        #Snap to nearest cell center
        dtmp <- sqrt((Scand[1] - grid[, 1])^2 + (Scand[2] - grid[, 2])^2)
        snap=which(dtmp==min(dtmp))
        Scand <- c(snap,as.numeric(grid[snap,1:2]))
        #business as usual
        dtmp <- sqrt((Scand[2] - X[, 1])^2 + (Scand[3] - X[, 2])^2)
        lamd.thisj<- lam0*exp(-D[j,]*D[j,]/(2*sigma*sigma))
        lamd.cand<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
        llS<- sum(func(lamd.thisj,y2D[j,],K,z[j],X))
        llcand<- sum(func(lamd.cand,y2D[j,],K,z[j],X))
        prior.s <- (beta0 + beta1 * (grid[s[j,1],3]==2))
        prior.s.cand <- (beta0 + beta1 * (grid[Scand[1],3]==2))
        if (runif(1) < exp((llcand+prior.s.cand) - (llS+prior.s))) {
          s[j, ] <- Scand
          D[j, ] <- dtmp
          lamd[j, ] <- lamd.cand
        }

      }

      #Do we record output on this iteration?
      if(i>nburn&i%%nthin==0){
        scellout[idx,]=s[,1]
        sxout[idx,]<- s[,2]
        syout[idx,]<- s[,3]
        zout[idx,]<- z
        out[idx,]<- c(lam0,sigma ,N,EN,beta0,beta1)
        idx=idx+1
      }
    }  # end of MCMC algorithm
    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout)
    }else{
      list(out=out)
    }
  }

#Fix upbeta
#
# out=upbeta(beta0, beta1, grid,  s[,1],proppars$beta0,proppars$beta1, EN,  N, cellArea, z)
# for(i in 1:100){
#   out=upbeta(out[[1]], out[[2]], grid, s[,1],proppars$beta0,proppars$beta1, out[[3]],  N, cellArea, z)
# }
