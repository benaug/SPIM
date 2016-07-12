N=50
p01=0.13
p02=0.26
lam01=-log(1-p01)
lam02=-log(1-p02)
1-exp(-c(lam01,lam02))
sigma=0.50
K=6
buff=2
X<- expand.grid(3:8,3:8)
X=cbind(X,2)
#Add number of detectors
t1=rep(1,nrow(X))
t1[which(X[,2]%in%c(3,4))]=2
t2=rep(1,nrow(X))
t2[which(X[,2]%in%c(5,6))]=2
t3=rep(1,nrow(X))
t3[which(X[,2]%in%c(7,8))]=2
tf=cbind(t1,t1,t2,t2,t3,t3)
#Simulate some data
data=sim2side(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff,tf=tf)
sim2side.plot(data,plottimes=c(1,3,5))
