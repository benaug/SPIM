library(SPIM)
N=50
p01=0.13
p02=0.26
lam01=-log(1-p01)
lam02=-log(1-p02)
1-exp(-c(lam01,lam02))
sigma=0.50
K=6
buff=2
niter=1000 #should run more than this and discard a burn in
nburn=1
nthin=1
xlim<- c(1,10)
ylim<- c(1,10)
X<- expand.grid(3:8,3:8)
#Add number of detectors
X=cbind(X,1)
X[which(X[,2]%in%c(4,7)),3]=2
#Simulate some data
data=sim2side(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff)
sim2side.plot(data,plottimes=c(1,3,5))

inits=list(psi=0.5,lam01=lam01,lam02=lam02,sigma=sigma)

#inits=list(psi1=0.5,psi2=0.5,psi3=0.5,lam01=lam01,lam02=lam02,sigma=sigma)



a=Sys.time()
store=mcmc.2side(data,niter=niter,nburn=nburn,nthin=nthin, M = 100,inits=inits,swap=10,keepACs=TRUE,Rcpp=TRUE)
#store=mcmc.2side.ind(data,niter=niter,nburn=nburn,nthin=nthin, M = 100,inits=inits,keepACs=TRUE,Rcpp=TRUE)
b=Sys.time()
b-a
p=pID(data,store$ID_L,store$ID_R)

#plot posteriors
par(mfrow=c(2,2))
plot(store$out[,1],main="lam01",type="l")
plot(store$out[,2],main="lam02",type="l")
plot(store$out[,3],main="sigma",type="l")
plot(store$out[,4],main="N",type="l")
par(mfrow=c(1,1))

#test build data
both2=data.frame(which(data$both[,1,,]>0,arr.ind=TRUE))
left2=data.frame(which(data$left[,2,,]>0,arr.ind=TRUE))
right2=data.frame(which(data$right[,3,,]>0,arr.ind=TRUE))
both2=data.frame(both2,type="B")
left2=data.frame(left2,type="L")
right2=data.frame(right2,type="R")
input=rbind(both2,left2,right2)
colnames(input)=c("ID","occ","trap","type")
data2=build.data(input,K,X,buff=2)
all(data2$both==data$both)
all(data2$left==data$left)
all(data2$right==data$right)
all(data2$IDknown==data$IDknown)


#test tf
data=sim2sidetf(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff,failprob=0.05,faildur=3)

psi=0.5
inits=list(psi=psi,lam01=lam01,lam02=lam02,sigma=sigma)
a=Sys.time()
store=mcmc.2side(data,niter=1000,nburn,nthin, inits,M = 120,Rcpp=TRUE)
b=Sys.time()
b-a



colMeans(store$out[500:999,])


##Test SCR0
N=100
p0=0.13
lam0=-log(1-p0)
sigma=0.50
K=5
buff=2
niter=1000 #should run more than this and discard a burn in
nburn=1
nthin=1
xlim<- c(1,10)
ylim<- c(1,10)
X<- expand.grid(3:8,3:8)

data=simSCR(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff)


inits=list(psi=0.5,lam0=lam0,sigma=sigma)
store=mcmc.SCR(data,niter=1000,nburn,nthin, inits,M = 120,Rcpp=TRUE)








library(coda)
effectiveSize(store$out)

##What to plot.  1.  estimated activity center of each real ID with it's real left and right captures.  Got the last part, need first part.
k=23 #left guy
i=data$ID_L[k] #Real individual
offset=0.1
plot(X[X[,3]==1,1:2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2,xlab="X",ylab="Y")
points(X[X[,3]==2,1]-offset,X[X[,3]==2,2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2)
points(X[X[,3]==2,1]+offset,X[X[,3]==2,2],xlim=xlim,ylim=ylim,pch=4,cex=1.5,lwd=2)
who=store$ID_Lout[,k]
for(j in 1:(niter-nburn)){
  points(store$sxout[j,who[j]],store$syout[j,who[j]],pch=20,col=rgb(1,0,0,0.1),cex=0.1)
}
points(data$s[data$IDknownR,],col="forestgreen",pch=20,cex=2)
points(data$s[(max(data$IDknownR)+1):N,],col="black",pch=20,cex=1)
singles=setdiff(unique(c(data$ID_L,data$ID_R)),1:nrow(data$Rboth))
points(data$s[singles,],col="gold",pch=20,cex=2)
points(data$s[i,1],data$s[i,2],col="red",pch=20,cex=2)

lefts=setdiff(data$ID_L,data$ID_R)
rights=setdiff(data$ID_R,data$ID_L)
LRs=setdiff(intersect(data$ID_R,data$ID_L),1:nrow(data$Rboth))
text(x=data$s[lefts,1],y=data$s[lefts,2],rep("L",length(lefts)),cex=0.75)
text(x=data$s[rights,1],y=data$s[rights,2],rep("R",length(lefts)),cex=0.75)
text(x=data$s[LRs,1],y=data$s[LRs,2],rep("LR",length(lefts)),cex=0.75)
text(x=data$s[data$IDknownR,1],y=data$s[data$IDknownR,2],rep("B",length(data$IDknownR)),cex=0.75)

#plot left caps
offset2=0.25
jit=0.1
l=which(data$ID_L==i)
if(length(l)>0){
  leftcaps=X[which(apply(data$Rleft,c(1,4),sum)[l,]>0),]
  text(x=jitter(leftcaps[,1],jit),y=jitter(leftcaps[,2]-offset2,jit),labels=rep("L",nrow(leftcaps)),cex=1.25,col="darkblue")
}
r=which(data$ID_R==i)
if(length(r)>0){
  rightcaps=X[which(apply(data$Rright,c(1,4),sum)[r,]>0),]
  text(x=jitter(rightcaps[,1],jit),y=jitter(rightcaps[,2]+offset2,jit),labels=rep("R",nrow(rightcaps)),cex=1.25,col="darkblue")
}
