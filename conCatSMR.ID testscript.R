source("sim.conCatSMR.ID.R")
source("mcmc.conCatSMR.ID.R")
library(coda)
N=50
n.marked=12
lam0=0.35
sigma=0.50
K=10 #number of occasions
buff=3 #state space buffer
X<- expand.grid(3:11,3:11) #make a trapping array
pMarkID=c(.8,.8)#probability of observing marked status of marked and unmarked individuals
pID=.8 #Probability marked individuals are identified
ncat=3  #number of ID categories
gamma=IDcobs=vector("list",ncat) #population frequencies of each category level. Assume equal here.
nlevels=rep(2,ncat) #number of levels per IDcat
for(i in 1:ncat){
  gamma[[i]]=rep(1/nlevels[i],nlevels[i])
  IDcovs[[i]]=1:nlevels[i]
}
pIDcat=rep(1,ncat)#category observation probabilities
tlocs=0 #telemetry locs/marked individual
obstype="poisson" #observation model, count or presence/absence?
marktype="premarked" #premarked or natural ID (marked individuals must be captured)?
data=sim.conCatSMR.ID(N=N,n.marked=n.marked,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,ncat=ncat,
                      pIDcat=pIDcat,gamma=gamma,IDcat=IDcat,pMarkID=pMarkID,tlocs=tlocs,pID=pID,marktype=marktype)
#MCMC
inits=list(lam0=lam0,sigma=sigma,gamma=gamma,psi=0.5) #start at simulated values
proppars=list(lam0=0.1,sigma=0.08,s=0.5,st=0.08) #st only for telemetered inds. Should be smaller than s.
M=150
keepACs=TRUE
keepGamma=FALSE
niter=1000
nburn=0
nthin=1
IDup="Gibbs"
out=mcmc.conCatSMR.ID(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,obstype=obstype,
           proppars=proppars,keepACs=TRUE,keepGamma=TRUE,IDup=IDup)


plot(mcmc(out$out))

1-rejectionRate(mcmc(out$out)) #shoot for 0.2 - 0.4 for lam0 and sigma. If too low, raise proppar. If too high, lower proppar.
1-rejectionRate(mcmc(out$sxout)) #activity center acceptance in x dimension. Shoot for min of 0.2
sum(rowSums(data$y.sight[(n.marked+1):N,,])>0) #true number of n.um
