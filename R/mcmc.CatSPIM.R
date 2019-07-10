#' Fit the categorical spatial partial identity model (catSPIM)
#' @description Will document after I defend. Should be able to get the idea in the example below.
#' @author Ben Augustine
#' @examples
#' \dontrun{library(coda)
#'#Normal SCR stuff
#'N=50
#'lam0=0.25
#'sigma=0.50
#'K=10
#'buff=3 #state space buffer. Should be at least 3 sigma.
#'X<- expand.grid(3:11,3:11)
#'obstype="poisson"
#'
#'#categorical identity covariate stuff
#'ncat=2  #number of ID covariates
#'gamma=vector("list",ncat) #population frequencies of each category level
#'
#'#Do this if you want to use the number of alleles at a locus to generate all possible genotypes
#'nlevels=rep(NA,ncat) #Number of levels per ID covariate
#'nallele=rep(3,ncat)  #number of alleles at each loci
#'for(i in 1:ncat){
#'  nlevels[i]=nallele[i]*(nallele[i]+1)/2
#'  gamma[[i]]=rep(1/nlevels[i],nlevels[i]) #generating all equal genotype frequencies
#'}
#'IDcovs=vector("list",ncat)#Store unique genotypes
#'for(i in 1:length(IDcovs)){
#'  IDcovs[[i]]=expand.grid(1:nallele[i],1:nallele[i])
#'  IDcovs[[i]]=unique(t(apply(IDcovs[[i]],1,sort)))
#'}
#'#sequentially number the unique genotypes
#'for(i in 1:ncat){
#'  IDcovs[[i]]=1:nrow(IDcovs[[i]])
#'}
#'
#'#Or this for generic ID covariates
#'gamma=vector("list",ncat)
#'nlevels=rep(3,ncat) #Number of levels per ID covariate
#'for(i in 1:ncat){
#'  gamma[[i]]=rep(1/nlevels[i],nlevels[i]) #generating all equal category level frequencies
#'}
#'IDcovs=vector("list",ncat)
#'for(i in 1:length(IDcovs)){
#'  IDcovs[[i]]=1:nlevels[i]
#'}
#'pID=rep(0.8,ncat)#sample by covariate level observation probability.  e.g. loci amplification probability
#'
#'#ncat=1 with nlevel=1 will produce unmarked SCR data with no ID covariates. 
#'#Well, everyone has the same covariate value so they are effectively unmarked
#'
#'#Simulate some data
#'data=simCatSPIM(N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,
#'ncat=ncat,pID=pID,gamma=gamma,
#'IDcovs=IDcovs)
#'str(data$y)#True data 
#'str(data$y.obs)#observed data
#'rowSums(data$y.obs)#one observation per row because they cannot be deterministically linked
#'str(data$G.true) #True categorical identities. Same # of rows as y
#'str(data$G.obs) #Observed categorical identities. Same # of rows as y.obs
#'data$G.obs #zeros are missing values
#'
#'#MCMC stuff
#'niter=2000 #how long to run the MCMC chain. 
#'nburn=0 #how much burnin to discard. I always do this afterwards so I can assess convergence better
#'nthin=1 #do we thin the chain. Only necessary if you need to reduce file size
#'nswap=nrow(data$y.obs)/2 #Number of latent IDs to swap on each iteration. Updating half seems to work fine.
#'#Theoretically a tradeoff between mixing and run time, but mixing is fine with not many swaps.
#'M=150 #Data augmentation level
#'inits=list(psi=0.3,lam0=lam0,sigma=sigma,gamma=gamma) #initial values. Using simulated values
#'proppars=list(lam0=0.05,sigma=0.075,sx=1,sy=1) #MCMC proposal parameters. Tuning these is the most difficult part.
#'#shoot for rejection rates between 0.2 and 0.4. If a parameter is not moving, the proposal parameter is too large.
#'IDup="MH" #Gibbs or metropolis-hastings latent ID update. Must use MH for binomial model.
#'#Both about equally as fast for Poisson
#'keepACs=TRUE #Do we store the activity center posteriors and other latent structure including the ID vector?
#'keepGamma=TRUE #Do we store the gamma posterior?
#'
#'a=Sys.time()
#'out=mcmc.CatSPIM(data,niter=niter,nburn=nburn,nthin=nthin, nswap=nswap,
#'                 M = M, inits=inits,proppars=proppars,obstype=obstype,
#'                 IDup=IDup,keepACs=keepACs,keepGamma=keepGamma)
#'b=Sys.time()
#'b-a
#'
#'#Let's see what happened!
#'plot(mcmc(out$out))
#'#The new, key parameter here is n, the number of individuals actually captured. 
#'#If it's not changing value, you have enough ID covariates to remove all uncertainty in n.
#'#Discard some of the burnin to see the posterior of n better
#'plot(mcmc(out$out[500:niter,]))
#'data$n #true value of n
#'
#'#If you kept the gamma posteriors, you can plot those, too.
#'plot(mcmc(out$gammaOut[[1]]))
#'plot(mcmc(out$gammaOut[[2]]))
#'gamma #True values
#'
#'#This will get you the acceptance probabilies. Can't change N, n, or psi.
#'1-rejectionRate(mcmc(out$out))
#'#This will get you the acceptance probabilities for the x dimension of the activity centers.
#'#Hard to tune because it will be lower when an activity center has samples allocated
#'#to it and higher when it does not. I think just make sure it's not under 0.1ish
#'1-rejectionRate(mcmc(out$sxout))
#'effectiveSize(out$out)#should shoot for at least 400 for N. 
#'
#'#OK, now try it with more/fewer ID covariates, category levels, or more missing data
#'
#'#Get posterior probability sample x and sample y came from same individual
#'niters=niter-nburn
#'
#'check=1 #Which sample to check
#'storematch=matrix(0,nrow=ncol(out$IDout),ncol=niters)
#'for(i in 1:niters){
#'  storematch[,i]=out$IDout[i,check]==out$IDout[i,]
#'}
#'rowSums(storematch)/niters
#'#True IDs stored here
#'data$IDtrue
#'}
#' @export

mcmc.CatSPIM <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=NA,obstype="poisson",nswap=NA,
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           keepACs=FALSE,keepGamma=FALSE,keepG=FALSE,IDup="Gibbs",priors=NA){
    ###
    library(abind)
    y.obs<-data$y.obs
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(y.obs)[3]
    ncat=data$IDlist$ncat
    nallele=data$IDlist$nallele
    IDcovs=data$IDlist$IDcovs
    buff<- data$buff
    xlim<- c(min(X[,1]),max(X[,1]))+c(-buff, buff)
    ylim<- c(min(X[,2]),max(X[,2]))+c(-buff, buff)
    n.samples=sum(y.obs)
    constraints=data$constraints
    G.obs=data$G.obs
    if(!is.matrix(G.obs)){
      G.obs=matrix(G.obs)
    }
    if(!is.list(IDcovs)){
      stop("IDcovs must be a list")
    }
    nlevels=unlist(lapply(IDcovs,length))
    if(ncol(G.obs)!=ncat){
      stop("G.obs needs ncat number of columns")
    }
    if(!all(is.na(priors))){
      if(!is.list(priors))stop("priors must be a list" )
      if(!all(names(priors)==c("sigma")))stop("priors list element name must be sigma")
      warning("Using Gamma prior for sigma")
      usePriors=TRUE
    }else{
      warning("No sigma prior entered, using uniform(0,infty).")
      usePriors=FALSE
    }

    
    #data checks
    if(length(dim(y.obs))!=3){
      stop("dim(y.obs) must be 3. Reduced to 2 during initialization")
    }

    if(is.na(nswap)){
      nswap=n.samples/2
      warning("nswap not specified, using n.samples/2")
    }
    if(!IDup%in%c("MH","Gibbs")){
      stop("IDup must be MH or Gibbs")
    }
    if(IDup=="MH"){
      # stop("MH not implemented, yet")
    }
    if(obstype=="bernoulli"&IDup=="Gibbs"){
      stop("Must use MH IDup for bernoulli data")
    }
    
    #If using polygon state space
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
    
    ##pull out initial values
    psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    gamma=inits$gamma
    if(!is.list(gamma)){
      stop("inits$gamma must be a list")
    }
    
    #make constraints if missing
    if(is.null(constraints)){
      constraints=matrix(1,nrow=n.samples,ncol=n.samples)
      for(i in 1:n.samples){
        for(j in 1:n.samples){
          guys1=which(G.obs[i,]!=0)
          guys2=which(G.obs[j,]!=0)
          comp=guys1[which(guys1%in%guys2)]
          if(any(G.obs[i,comp]!=G.obs[j,comp])){
            constraints[i,j]=0
          }
        }
      }
    }
    if(nrow(constraints)!=ncol(constraints)){
      stop("identity constraint matrix needs to be symmetric")
    }
    #If bernoulli data, add constraints that prevent y.true[i,j,k]>1
    binconstraints=FALSE
    if(obstype=="bernoulli"){
      # idx=which(y.obs>0,arr.ind=TRUE)
      idx=t(apply(y.obs,1,function(x){which(x>0,arr.ind=TRUE)}))
      for(i in 1:n.samples){
        for(j in 1:n.samples){
          if(i!=j){
            # if(all(idx[i,2:3]==idx[j,2:3])){
            if(all(idx[i,1:2]==idx[j,1:2])){
              constraints[i,j]=0 #can't combine samples from same trap and occasion in binomial model
              constraints[j,i]=0
              binconstraints=TRUE
            }
          }
        }
      }
    }
    
    #Build y.true
    y.true=array(0,dim=c(M,J,K))
    y.true2D=matrix(0,nrow=M,ncol=J)
    y.obs2D=apply(y.obs,c(1,2),sum)
    ID=rep(NA,n.samples)
    idx=1
    for(i in 1:n.samples){
      if(idx>M){
        stop("Need to raise M to initialize y.true")
      }
      traps=which(y.obs2D[i,]>0)
      # y.true2D=apply(y.true,c(1,2),sum)
      if(length(traps)==1){
        cand=which(y.true2D[,traps]>0)#guys caught at same traps
      }else{
        cand=which(rowSums(y.true2D[,traps])>0)#guys caught at same traps
      }
      if(length(cand)>0){
        if(length(cand)>1){#if more than 1 ID to match to, choose first one
          cand=cand[1]
        }
        #Check constraint matrix
        cands=which(ID%in%cand)#everyone assigned this ID
        if(all(constraints[i,cands]==1)){#focal consistent with all partials already assigned
          y.true[cand,,]=y.true[cand,,]+y.obs[i,,]
          ID[i]=cand
        }else{#focal not consistent
          y.true[idx,,]=y.obs[i,,]
          y.true2D[idx,traps]=y.true2D[idx,traps]+1
          ID[i]=idx
          idx=idx+1
        }
      }else{#no assigned samples at this trap
        y.true[idx,,]=y.obs[i,,]
        y.true2D[idx,traps]=y.true2D[idx,traps]+1
        ID[i]=idx
        idx=idx+1
      }
    }
    y.true2D=apply(y.true,c(1,2),sum)
    known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))
    if(binconstraints){
      if(any(y.true>1))stop("bernoulli data not initialized correctly")
    }
    #Check assignment consistency with constraints
    checkID=unique(ID)
    for(i in 1:length(checkID)){
      idx=which(ID==checkID[i])
      if(!all(constraints[idx,idx]==1)){
        stop("ID initialized improperly")
      }
    }
    y.true2D=apply(y.true,c(1,2),sum)
    known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))
    
    #Initialize z
    z=1*(apply(y.true2D,1,sum)>0)
    add=M*(0.5-sum(z)/M)
    if(add>0){
      z[sample(which(z==0),add)]=1 #switch some uncaptured z's to 1.
    }
    
    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.true)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y.true2D[i,]>0,1:2],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s[i,]<- trps
      }
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
    
    #collapse data to 2D
    y.true=y.true2D
    y.obs=apply(y.obs,c(1,2),sum)
    
    D=e2dist(s, X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    ll.y=array(0,dim=c(M,J))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      ll.y=dbinom(y.true,K,pd*z,log=TRUE)
    }else if(obstype=="poisson"){
      ll.y=dpois(y.true,K*lamd*z,log=TRUE)
    }
    ll.y.cand=ll.y
    
  
    
    #Initialize G.true
    G.true=matrix(0, nrow=M,ncol=ncat)
    for(i in 1:max(ID)){
      idx=which(ID==i)
      if(length(idx)==1){
        G.true[i,]=G.obs[idx,]
      }else{
        if(ncol(G.obs)>1){
          G.true[i,]=apply(G.obs[idx,],2, max) #consensus
        }else{
          G.true[i,]=max(G.obs[idx,])
        }
      }
    }
    G.latent=G.true==0
    for(j in 1:ncat){
      fix=G.true[,j]==0
      G.true[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
    }
    #unique IDcovs. 
    # genos=expand.grid(IDcovs)
    
    # some objects to hold the MCMC output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=5)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","n","psi"))
    if(keepACs){
      sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
      IDout=matrix(NA,nrow=nstore,ncol=length(ID))
    }
    idx=1 #for storing output not recorded every iteration
    if(keepGamma){
      gammaOut=vector("list",ncat)
      for(i in 1:ncat){
        gammaOut[[i]]=matrix(NA,nrow=nstore,ncol=nlevels[i])
        colnames(gammaOut[[i]])=paste("Lo",i,"G",1:nlevels[i],sep="")
      }
    }
    if(keepG){
      Gout=array(NA,dim=c(niter,M,ncat))
    }
    
    
    for(iter in 1:niter){
      #Update lam0
      if(obstype=="bernoulli"){
        llysum=sum(ll.y)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
          pd.cand=1-exp(-lamd.cand)
          ll.y.cand= dbinom(y.true,K,pd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
            llysum=llycandsum
          }
        }
        #Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
          pd.cand=1-exp(-lamd.cand)
          ll.y.cand= dbinom(y.true,K,pd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(usePriors){
            prior.curr=dgamma(sigma,priors$sigma[1],priors$sigma[2],log=TRUE)
            prior.cand=dgamma(sigma.cand,priors$sigma[1],priors$sigma[2],log=TRUE)
          }else{
            prior.curr=prior.cand=0
          }
          if(runif(1) < exp((llycandsum+prior.cand)-(llysum+prior.curr))){
            sigma<- sigma.cand
            lamd=lamd.cand
            pd=pd.cand
            ll.y=ll.y.cand
          }
        }
      }else{#poisson
        #Update lam0
        llysum=sum(ll.y)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
          ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
            llysum=llycandsum
          }
        }
        # Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
          ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(usePriors){
            prior.curr=dgamma(sigma,priors$sigma[1],priors$sigma[2],log=TRUE)
            prior.cand=dgamma(sigma.cand,priors$sigma[1],priors$sigma[2],log=TRUE)
          }else{
            prior.curr=prior.cand=0
          }
          if(runif(1) < exp((llycandsum+prior.cand)-(llysum+prior.curr))){
            sigma<- sigma.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
          }
        }
      }
      
      #ID update
      if(IDup=="Gibbs"){
        #Update y.true from full conditional canceling out inconsistent combos with constraints.
        up=sample(1:n.samples,nswap,replace=FALSE)
        for(l in up){
          nj=which(y.obs[l,]>0)
          #Can only swap if IDcovs match
          idx2=which(G.obs[l,]!=0)
          if(length(idx2)>1){#multiple loci observed
            possible=which(z==1&apply(G.true[,idx2],1,function(x){all(x==G.obs[l,idx2])}))
          }else if(length(idx2)==1){#single loci observed
            possible=which(z==1&G.true[,idx2]==G.obs[l,idx2])
          }else{#fully latent G.obs
            possible=which(z==1)#Can match anyone
          }
          njprobs=lamd[,nj]
          njprobs[setdiff(1:M,possible)]=0
          njprobs=njprobs/sum(njprobs)
          newID=sample(1:M,1,prob=njprobs)
          if(ID[l]!=newID){
            swapped=c(ID[l],newID)
            #update y.true
            y.true[ID[l],]=y.true[ID[l],]-y.obs[l,]
            y.true[newID,]=y.true[newID,]+y.obs[l,]
            ID[l]=newID
            ll.y[swapped,]= dpois(y.true[swapped,],K*lamd[swapped,],log=TRUE)
          }
        }
      }else{
        up=sample(1:n.samples,nswap,replace=FALSE)
        y.cand=y.true
        for(l in up){
          #find legal guys to swap with. z=1 and consistent constraints
          nj=which(y.obs[l,]>0)
          #Can only swap if IDcovs match
          idx2=which(G.obs[l,]!=0)
          if(length(idx2)>1){#multiple loci observed
            possible=which(z==1&apply(G.true[,idx2],1,function(x){all(x==G.obs[l,idx2])}))
          }else if(length(idx2)==1){#single loci observed
            possible=which(z==1&G.true[,idx2]==G.obs[l,idx2])
          }else{#fully latent G.obs
            possible=which(z==1)#Can match anyone
          }
          if(binconstraints){#can't have a y[i,j,k]>1
            legal=rep(TRUE,length(possible))
            for(i in 1:length(possible)){
              check=which(ID==possible[i])#Who else is currently assigned this possible new ID?
              if(length(check)>0){#if false, no samples assigned to this guy and legal stays true
                if(any(constraints[l,check]==0)){#if any members of the possible cluster are inconsistent with sample, illegal move
                  legal[i]=FALSE
                }
              }
            }
            possible=possible[legal]
          }
          njprobs=lamd[,nj]
          njprobs[setdiff(1:M,possible)]=0
          njprobs=njprobs/sum(njprobs)
          newID=ID
          newID[l]=sample(1:M,1,prob=njprobs)
          if(ID[l]==newID[l])next
          
          swapped=c(ID[l],newID[l])#order swap.out then swap.in
          propprob=njprobs[swapped[2]]
          backprob=njprobs[swapped[1]]
          #update y.true
          y.cand[ID[l],]=y.true[ID[l],]-y.obs[l,]
          y.cand[newID[l],]=y.true[newID[l],]+y.obs[l,]
          focalprob=(sum(ID==ID[l])/n.samples)*(y.true[ID[l],nj]/sum(y.true[ID[l],]))
          focalbackprob=(sum(newID==newID[l])/n.samples)*(y.cand[newID[l],nj]/sum(y.cand[newID[l],]))
          ##update ll.y
          if(obstype=="poisson"){
            ll.y.cand[swapped,]=dpois(y.cand[swapped,],K*lamd[swapped,],log=TRUE)
          }else{
            ll.y.cand[swapped,]=dbinom(y.cand[swapped,],K,pd[swapped,],log=TRUE)
          }
          if(runif(1)<exp(sum(ll.y.cand[swapped,])-sum(ll.y[swapped,]))*
             (backprob/propprob)*(focalbackprob/focalprob)){
            y.true[swapped,]=y.cand[swapped,]
            ll.y[swapped,]=ll.y.cand[swapped,]
            ID[l]=newID[l]
          }
        }
      }
      #update known.vector and G.latent
      known.vector=1*(rowSums(y.true)>0)
      G.true.tmp=matrix(0, nrow=M,ncol=ncat)
      for(i in unique(ID)){
        idx2=which(ID==i)
        if(length(idx2)==1){
          G.true.tmp[i,]=G.obs[idx2,]
        }else{
          if(ncol(G.obs)>1){
            G.true.tmp[i,]=apply(G.obs[idx2,],2, max) #consensus
          }else{
            G.true.tmp[i,]=max(G.obs[idx2,]) #consensus
          }
        }
      }
      G.latent=G.true.tmp==0
      
      #update G.true
      for(j in 1:ncat){
        swap=G.latent[,j]
        G.true[swap,j]=sample(IDcovs[[j]],sum(swap),replace=TRUE,prob=gamma[[j]])
      }
      
      #update genotype frequencies
      for(j in 1:ncat){
        x=rep(NA,nlevels[[j]])
        for(k in 1:nlevels[[j]]){
          x[k]=sum(G.true[z==1,j]==k)#genotype freqs in pop
        }
        gam=rgamma(rep(1,nlevels[[j]]),1+x)
        gamma[[j]]=gam/sum(gam)
      }
      
      # probability of not being captured in a trap AT ALL
      if(obstype=="poisson"){
        pd=1-exp(-lamd)
      }
      pbar=(1-pd)^K
      prob0<- exp(rowSums(log(pbar)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      if(obstype=="bernoulli"){
        ll.y= dbinom(y.true,K,pd*z,log=TRUE)
      }else{
        ll.y= dpois(y.true,K*lamd*z,log=TRUE)
      }
      
      
      psi=rbeta(1,1+sum(z),1+M-sum(z))
      
      ## Now we have to update the activity centers
      for (i in 1:M) {
        Scand <- c(rnorm(1, s[i, 1], proppars$sx), rnorm(1, s[i, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamd.cand[i,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          if(obstype=="bernoulli"){
            pd.cand[i,]=1-exp(-lamd.cand[i,])
            ll.y.cand[i,]= dbinom(y.true[i,],K,pd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              pd[i,]=pd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }else{#poisson
            ll.y.cand[i,]= dpois(y.true[i,],K*lamd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(keepACs){
          sxout[idx,]<- s[,1]
          syout[idx,]<- s[,2]
          zout[idx,]<- z
          IDout[idx,]=ID
        }
        if(keepGamma){
          for(k in 1:ncat){
            gammaOut[[k]][idx,]=gamma[[k]]
          }
        }
        if(keepG){
          Gout[idx,,]=G.true
        }
        out[idx,]<- c(lam0,sigma ,sum(z),length(unique(ID)),psi)
        idx=idx+1
      }
    }  # end of MCMC algorithm
    
    if(keepACs&keepGamma){
      out2=list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout,gammaOut=gammaOut)
    }else if(keepACs&!keepGamma){
      out2=list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else if(!keepACs&keepGamma){
      out2=list(out=out,gammaOut=gammaOut)
    }else{ 
      out2=list(out=out)
    }
    if(keepG){
      out2[[length(out2)+1]]=Gout
      names(out2)[length(out2)]="Gout"
    }
    out2
  }

