e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

cellprobsSCR<- function(lamd){
  # For gaussian hazard model convert lamda(s,x) to p(s,x)
  N<- dim(lamd)[1]
  J<- dim(lamd)[2] # traps
  # just 2 outcomes: left captures and right captures
  pmat<- matrix(NA,nrow=N,ncol=J)
  for(j in 1:J){
    pmat[,j]<- 1-exp(-lamd[,j])
  }
  pmat
}

simSCRhet <-
  function(N=120,lam0=0.2,sigma=0.70,sig.lam0=0,sig.sigma=0,sig.a0=0,beta=1,K=10,X=X,buff=3,compensate="Perfect",obstype="bernoulli"){
    #######Capture process######################
    # # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)

    #Heterogeneity stuff
    if(sig.lam0==0|sig.sigma==0){stop('must specify heterogeneity for lam0 or sigma')}
    if(compensate=="Perfect"){ #Perfect compensatory heterogeneity
      if(sig.a0==0){#no variation in a0
        sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
        a0=2*pi*sigma^2*lam0
        lam0.ind=a0/(2*pi*sigma.ind^2)
      }else{#variation in a0
        a0=2*pi*sigma^2*lam0
        a0.ind=exp(rnorm(N,log(a0),sig.a0))
        sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
        lam0.ind=a0.ind/(2*pi*sigma.ind^2)
      }
    }else if(compensate=="Imperfect"){#Deviation from perfect CH
      if(beta==1){stop('if compensate==Imperfect, beta must deviate from 1')}
      if(sig.a0==0){#no variation in a0
        sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
        a0=2*pi*sigma^2*lam0*beta
        lam0.ind=a0/(2*pi*sigma.ind^2)
      }else{#variation in a0
        a0=2*pi*sigma^2*lam0*beta
        a0.ind=exp(rnorm(N,log(a0),sig.a0))
        sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
        lam0.ind=a0.ind/(2*pi*sigma.ind^2)
      }

    }else if(compensate=="Independent"){
      if(sig.lam0==0&sig.sigma==0){stop('must specify heterogeneity for lam0 and sigma if compansate==Independent')}
      sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
      lam0.ind=exp(rnorm(N,log(lam0),sig.lam0))
    }else{
      stop("compensate must equal Perfect or Imperfect")

    }

    beta=0.25
    a0=2*pi*sigma^2*lam0*beta
    sig.a0=0.2
    a0.ind=exp(rnorm(N,log(a0),sig.a0))
    sigma.ind=exp(rnorm(N,log(sigma),sig.sigma))
    lam.ind=a0.ind/(2*pi*sigma.ind^2)
    plot(lam.ind~sigma.ind)
    plot()


    # Simulate encounter history
    y <-array(0,dim=c(N,K,J))
    if(obstype=="bernoulli"){
      pd=cellprobsSCR(lamd)
      if(lam0b==0){ #if no behavioral response
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,pd[i,j])
            }
          }
        }
      }else{
        if((-lam0b)>lam0){stop("-b must be > lam0")}
        lamdb=(lam0+lam0b)*exp(-D*D/(2*sigma*sigma))
        pdb=cellprobsSCR(lamdb)
        state=matrix(0,nrow=N,ncol=J) #Matrix of indices 1 indicating previously captured at trap 0 o.w.
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,pd[i,j]*(1-state[i,j])+pdb[i,j]*state[i,j])
              state[i,j]=max(state[i,j],y[i,k,j])  #update state
            }
          }
        }
      }
    }else if(obstype=="poisson"){
      if(lam0b==0){ #if no behavioral response
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rbinom(1,1,lamd[i,j])
            }
          }
        }
      }else{
        if((-lam0b)>lam0){stop("-b must be > lam0")}
        lamdb=(lam0+lam0b)*exp(-D*D/(2*sigma*sigma))
        state=matrix(0,nrow=N,ncol=J) #Matrix of indices 1 indicating previously captured at trap 0 o.w.
        for(i in 1:N){
          for(j in 1:J){
            for(k in 1:K){
              y[i,k,j]=rpois(1,lamd[i,j]*(1-state[i,j])+lamdb[i,j]*state[i,j])
              state[i,j]=max(state[i,j],y[i,k,j])  #update cap
            }
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }
    caps=apply(y,1,sum)
    idx=order(caps,decreasing=TRUE)
    y=y[idx,,]
    s=s[idx,]
    n=sum(caps>0)
    out<-list(y=y,s=s,X=X, K=K,n=n,buff=buff,obstype=obstype)
    return(out)
  }
