library(SPIM)
library(coda)
source("simSCR2DNA.R")
source("SCR2DNAmcmc.R")
N=50
p01=0.2
p02=0.2
lam01=-log(1-p01)
lam02=-log(1-p02)
1-exp(-c(lam01,lam02))
sigma=0.50
K1=12
K2=10
#Hair snare grid
xlim<- c(1,10)
ylim<- c(1,10)
# X1<- expand.grid(3:8,3:8)
X1=cbind(3:8,8)
X1=rbind(X1,cbind(3:8,3))
X1=rbind(X1,cbind(3:8,5))
X1=rbind(X1,cbind(3:8,6))

#Scat grid
X2<-cbind(seq(3.5,7.5,0.5),3.5)
X2=rbind(X2,cbind(seq(3.5,7.5,0.5),4.5))
X2=rbind(X2,cbind(seq(3.5,7.5,0.5),5.5))
X2=rbind(X2,cbind(seq(3.5,7.5,0.5),6.5))
X2=rbind(X2,cbind(seq(3.5,7.5,0.5),7.5))
buff=2


data=simSCR2DNA(N=N,lam01=lam01,lam01b=0,lam02=lam02,sigma=sigma,K1=K1,K2=K2,X1=X1,X2=X2,buff=buff)
inits=list(psi=0.5,lam01=lam01,lam02=lam02,sigma=sigma)
proppars=list(lam01=0.1,lam02=0.1,sigma=0.05,sx=0.2,sy=0.2)
niter=5000
nburn=1000
nthin=1
M=80
store=SCR2DNAmcmc(data,niter=niter,nburn=nburn, nthin=nthin, M =M, inits=inits,proppars=proppars,keepACs=TRUE)

plot(mcmc(store$out))
summary(mcmc(store$out))


both=which(rowSums(data$y1)>0&rowSums(data$y2)>0)
hair=which(rowSums(data$y1)>0&rowSums(data$y2)==0)
scat=which(rowSums(data$y1)==0&rowSums(data$y2)>0)
none=which(rowSums(data$y1)==0&rowSums(data$y2)==0)

plot(X1,xlim=xlim,ylim=ylim,pch=4)
points(X2)
points(data$s[both,],pch=16,col="green")
points(data$s[hair,],pch=16,col="red")
points(data$s[scat,],pch=16,col="blue")
points(data$s[none,],pch=16,col="black")

