#' Run MCMC algorithm for basic SCR model with a trap functionality file using Rcpp.
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
SCRmcmctfRcpp <-
function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=FALSE){
###
library(abind)
y<-data$y
X<-as.matrix(data$X)
J<-nrow(X)
K<- dim(y)[2]
n<- data$n


##pull out initial values
psi<- inits$psi
lam0<- inits$lam0
sigma<- inits$sigma
if("vertices"%in%names(data)){
  vertices=data$vertices
  useverts=TRUE
  xlim=c(min(vertices[,1]),max(vertices[,1]))
  ylim=c(min(vertices[,2]),max(vertices[,2]))
}else{
  useverts=FALSE
  buff<- data$buff
  xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
  ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
  vertices=rbind(xlim,ylim)
}
#augment data
y<- abind(y,array(0, dim=c( M-dim(y)[1],K, J)), along=1)
known.vector=c(rep(1,data$n),rep(0,M-data$n))

#Make initial complete data set
y2D=apply(y,c(1,3),sum)
z=1*(apply(y2D,1,sum)>0)
z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  half is arbitrary. smarter way?

#trap history
tf=data$tf
if(length(dim(tf))==2){
  Ktf=rowSums(tf)
}else{
  Ktf=tf
}
K2D=matrix(rep(Ktf,M),nrow=M,ncol=J,byrow=TRUE)

#Optimize starting locations given where they are trapped.
s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
idx=which(rowSums(y)>0) #switch for those actually caught
for(i in idx){
  trps<- X[y2D[i,]>0,1:2]
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


# some objects to hold the MCMC simulation output
nstore=(niter-nburn)/nthin
if(nburn%%nthin!=0){
  nstore=nstore+1
}
out<-matrix(NA,nrow=nstore,ncol=3)
dimnames(out)<-list(NULL,c("lam0","sigma","N"))
sxout<- syout<- zout<- ID_Lout<- ID_Rout<-matrix(NA,nrow=nstore,ncol=M)
idx=1 #for storing output not recorded every iteration

D<- e2dist(s, X)
lamd<- lam0*exp(-D*D/(2*sigma*sigma))

#Run MCMC
store=SCRRcpp::MCMC1tf( lam0,  sigma,y2D, z,  X, Ktf,D,known.vector,s,psi,xlim,ylim,useverts,vertices,proppars$lam0,proppars$sigma,proppars$sx,
                        proppars$sy,niter,nburn,nthin)
if(keepACs){
  list(out=store[[1]], sxout=store[[2]], syout=store[[3]], zout=store[[4]])
}else{
  list(out=store[[1]])
}
}
