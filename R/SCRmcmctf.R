#' Run MCMC algorithm for basic SCR model with a trap functionality file.
#' @param data a list produced by simSCR or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam0=lam0,sigma=sigma)
#' @param proppars a list of tuning parameters for the proposal distributions
#' @return  a list with the posteriors for the SCR parameters (out), s, z
#' @author Ben Augustine, Andy Royle
#' @description This function runs the MCMC algorithm for the basic SCR model.  The data list should have the following elements:
#' 1.  y, a n x J capture history
#' 2.  X,  a matrix with the X and Y trap locations in the first two columns and the number of cameras (1 or 2) at each trap in the third.
#' 3. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' 4.  tf, a vector of length J containing the number of occasions each trap was operational
#' @export

SCRmcmctf <-
function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2)){
###
library(abind)
y<-data$y
X<-as.matrix(data$X)
J<-nrow(X)
K<- dim(y)[2]
n<- dim(y)[1]

tf=data$tf
##pull out initial values
psi<- inits$psi
lam0<- inits$lam0
sigma<- inits$sigma

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
#augment data
y<- abind(y,array(0, dim=c( M-dim(y)[1],K, J)), along=1)
known.vector=c(rep(1,data$n),rep(0,M-data$n))

#trap history
if(length(dim(tf))==2){
  Ktf=rowSums(tf)
}else{
  Ktf=tf
}
K2D=matrix(rep(Ktf,M),nrow=M,ncol=J,byrow=TRUE)

#Make initial complete data set
y2D=apply(y,c(1,3),sum)
z=1*(apply(y2D,1,sum)>0)
z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?

#Optimize starting locations given where they are trapped.
s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx=which(rowSums(y)>0) #switch for those actually caught
for(i in idx){
  trps<- X[y2D[i,]>0,1:2]
  trps<-matrix(trps,ncol=2,byrow=FALSE)
  s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
}
#check to make sure everyone is in polygon
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


#Bernoulli Likelihood function
func<- function(lamd,y,K2D,z,X){
  #convert lamd to pd (gaussian hazard model)
  pd=1-exp(-lamd)
  #If data is M x K
  if(is.matrix(y)){
      v <-  dbinom(y,K2D,pd,log=TRUE)
      v[z==0,]<- 0
  }else{
  #If data is 1 x K
    v <- dbinom(y,K2D,pd,log=TRUE)
    v<- v*z
  }
  v
}

# some objects to hold the MCMC simulation output
nstore=(niter-nburn)/nthin
if(nburn%%nthin!=0){
  nstore=nstore+1
}
out<-matrix(NA,nrow=nstore,ncol=3)
dimnames(out)<-list(NULL,c("lam0","sigma","N"))
sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
idx=1 #for storing output not recorded every iteration

D<- e2dist(s, X)
lamd<- lam0*exp(-D*D/(2*sigma*sigma))

for(i in 1:niter){
  #Update lam0
  lik.curr<-  sum( func(lamd,y2D,K2D,z,X) )
    lam0.cand<- rnorm(1,lam0,proppars$lam0)
    if(lam0.cand > 0){
      lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
      lik.new<-  sum( func(lamd.cand,y2D,K2D,z,X) )
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
    lik.new<-  sum( func(lamd.cand,y2D,K2D,z,X) )
    if(runif(1) < exp(lik.new -lik.curr)){
      sigma<- sigma.cand
      lamd=lamd.cand
      lik.curr<- lik.new
    }
  }
  #Update psi gibbs
  ## probability of not being captured in a trap AT ALL
  pd=1-exp(-lamd)
  pbar=(1-pd)^K2D
  prob0<- exp(rowSums(log(pbar)))
  fc<- prob0*psi/(prob0*psi + 1-psi)
  z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
  lik.curr<-  sum( func(lamd,y,K2D,z,X) )
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
      lamd.thisj<- lam0*exp(-D[j,]*D[j,]/(2*sigma*sigma))
      lamd.cand<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
      llS<- sum(func(lamd.thisj,y2D[j,],Ktf,z[j],X))
      llcand<- sum(func(lamd.cand,y2D[j,],Ktf,z[j],X))
      if (runif(1) < exp(llcand - llS)) {
        s[j, ] <- Scand
        D[j, ] <- dtmp
        lamd[j, ] <- lamd.cand
      }
    }
  }
  #Do we record output on this iteration?
  if(i>nburn&i%%nthin==0){
    sxout[idx,]<- s[,1]
    syout[idx,]<- s[,2]
    zout[idx,]<- z
    out[idx,]<- c(lam0,sigma ,sum(z))
    idx=idx+1
  }
}  # end of MCMC algorithm

list(out=out, sxout=sxout, syout=syout, zout=zout)
}

