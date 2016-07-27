SCRmcmcRcpp <-
function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
###
library(abind)
y<-data$y
X<-as.matrix(data$X)
J<-nrow(X)
n<- dim(y)[1]

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

#Augment data and make initial complete data set
if(length(dim(y))==3){
  idx=which(rowSums(y)==0)
  if(length(idx)>0){
    y=y[-idx,,]
  }
  K<- dim(y)[2]
  y<- abind(y,array(0, dim=c( M-dim(y)[1],K, J)), along=1)
  y2D=apply(y,c(1,3),sum)
}else if(length(dim(y)==2)){
  if(is.na(K)){
    stop("if y is 2D, must supply K")
  }
  if(length(idx)>0){
    idx=which(rowSums(y)==0)
    y=y[-idx,]
  }
  y2D=y<- abind(y,array(0, dim=c( M-dim(y)[1],J)), along=1)
}else{
  stop("y must be either 2D or 3D")
}
known.vector=c(rep(1,data$n),rep(0,M-data$n))
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
store=SPIM::MCMC1( lam0,  sigma,y2D, z,  X, K,D,known.vector,s,psi,xlim,ylim,useverts,vertices,proppars$lam0,proppars$sigma,proppars$sx,
                      proppars$sy,niter,nburn,nthin)

if(keepACs){
  list(out=store[[1]], sxout=store[[2]], syout=store[[3]], zout=store[[4]],data=data)
}else{
  list(out=store[[1]],data=data)
}
}
