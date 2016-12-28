SCR2DNAMbmcmc <-
function(data,niter=2400,nburn=1200, nthin=5, M = 200, sharesig=FALSE,inits=inits,proppars=list(lam01=0.05,lam01b=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
###
  if(sharesig==FALSE){
    if(length(proppars$sigma)!=2|length(inits$sigma)!=2){
      stop("must supply 2 starting values and proppars if sharesig=FALSE")
    }
  }else{
    if(length(proppars$sigma)!=1|length(inits$sigma)!=1){
      stop("must supply only 1 starting value and proppars if sharesig=TRUE")
    }
    inits$sigma=rep(inits$sigma,2)
  }

  library(abind)
  y1<-data$y1
  y2<-data$y2
  X1<-as.matrix(data$X1)
X2<-as.matrix(data$X2)
J1<-nrow(X1)
J2<-nrow(X2)
#Remove guys not captured.
rem=which(rowSums(y1)==0&rowSums(y2)==0)
y1=y1[-rem,,]
y2=y2[-rem,,]
n<- dim(y1)[1]

#If using polygon state space
if("vertices"%in%names(data)){
  vertices=data$vertices
  useverts=TRUE
  xlim=c(min(vertices[,1]),max(vertices[,1]))
  ylim=c(min(vertices[,2]),max(vertices[,2]))
}else if("buff"%in%names(data)){
  buff<- data$buff
  xlim<- c(min(c(X1[,1],X2[,1])),max(c(X1[,1],X2[,1])))+c(-buff, buff)
  ylim<- c(min(c(X1[,2],X2[,2])),max(c(X1[,2],X2[,2])))+c(-buff, buff)
  vertices=cbind(xlim,ylim)
  useverts=FALSE
}else{
  stop("user must supply either 'buff' or 'vertices' in data object")
}
##pull out initial values
psi<- inits$psi
lam01<- inits$lam01
lam01b<- inits$lam01b
lam02<- inits$lam02
sigma<- inits$sigma

#Augment data and make initial complete data set
#Data set 1
if(length(dim(y1))==3){
  K1<- dim(y1)[2]
  y1<- abind(y1,array(0, dim=c( M-dim(y1)[1],K1, J1)), along=1)
  y12D=apply(y1,c(1,3),sum)
}else if(length(dim(y1)==2)){
  if(is.na(K)){
    stop("if y is 2D, must supply K")
  }
  y12D=abind(y1,array(0, dim=c( M-dim(y1)[1],J1)), along=1)
}else{
  stop("y must be either 2D or 3D")
}
#Data set 2
if(length(dim(y2))==3){
  K2<- dim(y2)[2]
  y2<- abind(y2,array(0, dim=c( M-dim(y2)[1],K2, J2)), along=1)
  y22D=apply(y2,c(1,3),sum)
}else if(length(dim(y2)==2)){
  if(is.na(K)){
    stop("if y is 2D, must supply K")
  }
  y22D=abind(y2,array(0, dim=c( M-dim(y2)[1],J2)), along=1)
}else{
  stop("y must be either 2D or 3D")
}


known.vector=c(rep(1,n),rep(0,M-n))

z=known.vector
z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?

#Optimize starting locations given where they are trapped.
s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx=which(known.vector==1) #switch for those actually caught
for(i in idx){
  trps<- c(X1[y12D[i,]>0,1:2],X2[y22D[i,]>0,1:2])
  trps<-matrix(trps,ncol=2,byrow=FALSE)
  s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
}
#check to make sure everyone is in polygon
if("vertices"%in%names(data)){
  vertices=data$vertices
  useverts=TRUE
}else{
  useverts=FALSE
}
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

#Behavioral state matrix
state=array(0,dim=dim(y1))
for(i in 1:N){
  for(k in 1:(K1-1)){
    state[i,k+1,which(y1[i,k,]==1|state[i,k,]==1)]=1
  }
}

#Bernoulli Likelihood function
func<- function(lamd1,lamd1b,lamd2,y1,y2,state,K2,z){
  #convert lamd to pd (gaussian hazard model)
  pd1=1-exp(-lamd1)
  pd1b=1-exp(-lamd1b)
  pd2=1-exp(-lamd2)
  #If data is M x K
  if(length(dim(y1))==3){
    v <-  apply(dbinom(y1,1,pd1*(state==0)+pd1b*(state==1),log=TRUE),1,sum)+rowSums(dbinom(y2,K2,pd2,log=TRUE))
    v[z==0]<- 0
  }else{
  #If data is 1 x K
    v <- sum(dbinom(y1,1,pd1*(state==0)+pd1b*(state==1),log=TRUE))+sum(dbinom(y2,K2,pd2,log=TRUE))
    v<- v*z
  }
  v
}


# some objects to hold the MCMC simulation output
nstore=(niter-nburn)/nthin
if(nburn%%nthin!=0){
  nstore=nstore+1
}
if(sharesig==TRUE){
  out<-matrix(NA,nrow=nstore,ncol=5)
  dimnames(out)<-list(NULL,c("lam01","lam01b","lam02","sigma","N"))
}else{
  out<-matrix(NA,nrow=nstore,ncol=6)
  dimnames(out)<-list(NULL,c("lam01","lam01b","lam02","sigma1","sigma2","N"))
}
sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
idx=1 #for storing output not recorded every iteration

D1<- e2dist(s, X1)
D2<- e2dist(s, X2)
lamd1<- lam01*exp(-D1*D1/(2*sigma[1]*sigma[1]))
lamd1=do.call(abind, c(list(replicate(K1,lamd1))))
lamd1=aperm(lamd1,c(1,3,2))
lamd1b<- lam01b*exp(-D1*D1/(2*sigma[1]*sigma[1]))
lamd1b=do.call(abind, c(list(replicate(K1,lamd1b))))
lamd1b=aperm(lamd1b,c(1,3,2))
lamd2<- lam02*exp(-D2*D2/(2*sigma[2]*sigma[2]))

for(i in 1:niter){
  #Update lam01
  lik.curr<-  sum( func(lamd1,lamd1b,lamd2,y1,y22D,state,K2,z) )
  lam01.cand<- rnorm(1,lam01,proppars$lam01)
  if(lam01.cand > 0){
    lamd1.cand<- lam01.cand*exp(-D1*D1/(2*sigma[1]*sigma[1]))
    lamd1.cand=do.call(abind, c(list(replicate(K1,lamd1.cand))))
    lamd1.cand=aperm(lamd1.cand,c(1,3,2))
    lik.new<-  sum( func(lamd1.cand,lamd1b,lamd2,y1,y22D,state,K2,z) )
    if(runif(1) < exp(lik.new -lik.curr)){
      lam01<- lam01.cand
      lamd1=lamd1.cand
      lik.curr<- lik.new
    }
  }
  #Update lam01b
  lam01b.cand<- rnorm(1,lam01b,proppars$lam01b)
  if(lam01b.cand > 0){
    lamd1b.cand<- lam01b.cand*exp(-D1*D1/(2*sigma[1]*sigma[1]))
    lamd1b.cand=do.call(abind, c(list(replicate(K1,lamd1b.cand))))
    lamd1b.cand=aperm(lamd1b.cand,c(1,3,2))
    lik.new<-  sum( func(lamd1,lamd1b.cand,lamd2,y1,y22D,state,K2,z) )
    if(runif(1) < exp(lik.new -lik.curr)){
      lam01b<- lam01b.cand
      lamd1b=lamd1b.cand
      lik.curr<- lik.new
    }
  }
  #Update lam02
  lam02.cand<- rnorm(1,lam02,proppars$lam02)
  if(lam02.cand > 0){
    lamd2.cand<- lam02.cand*exp(-D2*D2/(2*sigma[2]*sigma[2]))
    lik.new<-  sum( func(lamd1,lamd1b,lamd2.cand,y1,y22D,state,K2,z) )
    if(runif(1) < exp(lik.new -lik.curr)){
      lam02<- lam02.cand
      lamd2=lamd2.cand
      lik.curr<- lik.new
    }
  }
  #Update sigma
  if(sharesig==TRUE){
    sigma.cand<- rnorm(1,sigma,proppars$sigma)
    if(sigma.cand > 0){
      lamd1.cand<- lam01*exp(-D1*D1/(2*sigma.cand*sigma.cand))
      lamd1.cand=do.call(abind, c(list(replicate(K1,lamd1.cand))))
      lamd1.cand=aperm(lamd1.cand,c(1,3,2))
      lamd1b.cand<- lam01b*exp(-D1*D1/(2*sigma.cand*sigma.cand))
      lamd1b.cand=do.call(abind, c(list(replicate(K1,lamd1b.cand))))
      lamd1b.cand=aperm(lamd1b.cand,c(1,3,2))
      lamd2.cand<- lam02*exp(-D2*D2/(2*sigma.cand*sigma.cand))
      lik.new<-   sum( func(lamd1.cand,lamd1b.cand,lamd2.cand,y1,y22D,state,K2,z) )
      if(runif(1) < exp(lik.new -lik.curr)){
        sigma<- rep(sigma.cand,2)
        lamd1=lamd1.cand
        lamd1b=lamd1b.cand
        lamd2=lamd2.cand
        lik.curr<- lik.new
      }
    }
  }else{
    #update sigma 1
    sigma.cand<- rnorm(1,sigma[1],proppars$sigma[1])
    if(sigma.cand > 0){
      lamd1.cand<- lam01*exp(-D1*D1/(2*sigma.cand*sigma.cand))
      lamd1.cand=do.call(abind, c(list(replicate(K1,lamd1.cand))))
      lamd1.cand=aperm(lamd1.cand,c(1,3,2))
      lamd1b.cand<- lam01b*exp(-D1*D1/(2*sigma.cand*sigma.cand))
      lamd1b.cand=do.call(abind, c(list(replicate(K1,lamd1b.cand))))
      lamd1b.cand=aperm(lamd1b.cand,c(1,3,2))
      lik.new<-   sum( func(lamd1.cand,lamd1b.cand,lamd2,y1,y22D,state,K2,z) )
      if(runif(1) < exp(lik.new -lik.curr)){
        sigma[1]<- sigma.cand
        lamd1=lamd1.cand
        lamd1b=lamd1b.cand
        lik.curr<- lik.new
      }
    }
    #update sigma 2
    sigma.cand<- rnorm(1,sigma[2],proppars$sigma[2])
    if(sigma.cand > 0){
      lamd2.cand<- lam02*exp(-D2*D2/(2*sigma.cand*sigma.cand))
      lik.new<-   sum( func(lamd1,lamd1b,lamd2.cand,y1,y22D,state,K2,z) )
      if(runif(1) < exp(lik.new -lik.curr)){
        sigma[2]<- sigma.cand
        lamd2=lamd2.cand
        lik.curr<- lik.new
      }
    }
  }
  #Update psi gibbs
  ## probability of not being captured in a trap AT ALL
  pd1=1-exp(-lamd1[,1,])
  pd2=1-exp(-lamd2)
  pbar1=(1-pd1)^K1
  pbar2=(1-pd2)^K2
  prob0<- exp(rowSums(log(pbar1))+rowSums(log(pbar2)))
  fc<- prob0*psi/(prob0*psi + 1-psi)
  z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
  lik.curr<-   sum( func(lamd1,lamd1b,lamd2,y1,y22D,state,K2,z) )
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
      d1tmp <- sqrt((Scand[1] - X1[, 1])^2 + (Scand[2] - X1[, 2])^2)
      lamd1.thisj<- lamd1[j,,]
      lamd1.cand<- lam01*exp(-d1tmp*d1tmp/(2*sigma[1]*sigma[1]))
      lamd1.cand=t(do.call(rbind, c(list(replicate(K1,lamd1.cand)))))
      lamd1b.thisj<- lamd1b[j,,]
      lamd1b.cand<- lam01b*exp(-d1tmp*d1tmp/(2*sigma[1]*sigma[1]))
      lamd1b.cand=t(do.call(rbind, c(list(replicate(K1,lamd1b.cand)))))
      d2tmp <- sqrt((Scand[1] - X2[, 1])^2 + (Scand[2] - X2[, 2])^2)
      lamd2.thisj<- lam02*exp(-D2[j,]*D2[j,]/(2*sigma[2]*sigma[2]))
      lamd2.cand<- lam02*exp(-d2tmp*d2tmp/(2*sigma[2]*sigma[2]))
      llS<- sum(func(lamd1.thisj,lamd1b.thisj,lamd2.thisj,y1[j,,],y22D[j,],state[j,,],K2,z[j]))
      llcand<- sum(func(lamd1.cand,lamd1b.cand,lamd2.cand,y1[j,,],y22D[j,],state[j,,],K2,z[j]))
      if (runif(1) < exp(llcand - llS)) {
        s[j, ] <- Scand
        D1[j, ] <- d1tmp
        D2[j, ] <- d2tmp
        lamd1[j,, ] <- lamd1.cand
        lamd1b[j,, ] <- lamd1b.cand
        lamd2[j, ] <- lamd2.cand
      }
    }
  }
  #Do we record output on this iteration?
  if(i>nburn&i%%nthin==0){
    sxout[idx,]<- s[,1]
    syout[idx,]<- s[,2]
    zout[idx,]<- z
    if(sharesig==TRUE){
      out[idx,]<- c(lam01,lam01b,lam02,sigma[1],sum(z))
    }else{
      out[idx,]<- c(lam01,lam01b,lam02,sigma,sum(z))
    }
    idx=idx+1
  }
}  # end of MCMC algorithm

if(keepACs==TRUE){
  list(out=out, sxout=sxout, syout=syout, zout=zout)
}else{
  list(out=out)
}
}

