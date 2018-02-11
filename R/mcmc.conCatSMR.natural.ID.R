#' Fit the categorical spatial mark resight model for natural marks or other situation where number of marked inds is unknown
#' @description Will document after I defend. Should be able to get the idea in the example below.
#' @author Ben Augustine
#' @examples
#' \dontrun{library(coda)
#' N=50
#' n.marked=12
#' lam0=0.35
#' sigma=0.50
#' K=10 #number of occasions
#' buff=3 #state space buffer
#' X<- expand.grid(3:11,3:11) #make a trapping array
#' pMarkID=c(.8,.8)#probability of observing marked status of marked and unmarked individuals
#' pID=.8 #Probability marked individuals are identified
#' ncat=3  #number of ID categories
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each category level. Assume equal here.
#' nlevels=rep(2,ncat) #number of levels per IDcat
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#category observation probabilities
#' tlocs=0 #telemetry locs/marked individual
#' obstype="poisson" #observation model, count or presence/absence?
#' marktype="natural" #premarked or natural ID (marked individuals must be captured)?
#' data=sim.conCatSMR.ID(N=N,n.marked=n.marked,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,ncat=ncat,
#'                       pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs,pID=pID,marktype=marktype)
#' #MCMC
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,psi1=0.5,psi2=0.5) #start at simulated values
#' proppars=list(lam0=0.1,sigma=0.08,s=0.5,st=0.08) #st only for telemetered inds. Should be smaller than s.
#' M1=50 #marked data augmentation
#' M2=100 #unmarked data augmentation
#' keepACs=TRUE
#' keepGamma=FALSE
#' niter=1000
#' nburn=0
#' nthin=1
#' IDup="Gibbs"
#' out=mcmc.conCatSMR.natural.ID(data,niter=niter,nburn=nburn, nthin=nthin, M1 = M1,M2=M2, inits=inits,obstype=obstype,
#'                               proppars=proppars,keepACs=TRUE,keepGamma=TRUE,IDup=IDup)
#' 
#' 
#' plot(mcmc(out$out))
#' 
#' 1-rejectionRate(mcmc(out$out)) #shoot for 0.2 - 0.4 for lam0 and sigma. If too low, raise proppar. If too high, lower proppar.
#' 1-rejectionRate(mcmc(out$s1xout)) #marked activity center acceptance in x dimension. Shoot for min of 0.2
#' 1-rejectionRate(mcmc(out$s2xout)) #unmarked activity center acceptance in x dimension. Shoot for min of 0.2}
#' @export
mcmc.conCatSMR.natural.ID <-
  function(data,niter=2400,nburn=1200, nthin=5, M1 = 30,M2=200, inits=NA,obstype="poisson",nswap=NA,
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           keepACs=TRUE,keepGamma=TRUE,IDup="Gibbs"){
    ###
    library(abind)
    y.sight.marked=data$y.sight.marked
    y.sight.unmarked=data$y.sight.unmarked
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(y.sight.marked)[3]
    ncat=data$IDlist$ncat
    IDcovs=data$IDlist$IDcovs
    buff<- data$buff
    marked.guys=1:dim(y.sight.marked)[1]
    G.marked=data$G.marked
    G.unmarked=data$G.unmarked
  
    if(!is.matrix(G.marked)){
      G.marked=matrix(G.marked)
    }
    if(!is.matrix(G.unmarked)){
      G.unmarked=matrix(G.unmarked)
    }
    if(!is.list(IDcovs)){
      stop("IDcovs must be a list")
    }
    nIDcovs=unlist(lapply(IDcovs,length))
    if(ncol(G.marked)!=ncat){
      stop("G.marked needs ncat number of columns")
    }
    if(ncol(G.unmarked)!=ncat){
      stop("G.unmarked needs ncat number of columns")
    }
    #Are there unknown marked status guys?
    useUnk=FALSE
    if("G.unk"%in%names(data)){
      if(!is.na(data$G.unk[1])){
        G.unk=data$G.unk
        if(!is.matrix(G.unk)){
          G.marked=matrix(G.unk)
        }
        if(ncol(G.unk)!=ncat){
          stop("G.unk needs ncat number of columns")
        }
        y.sight.unk=data$y.sight.unk
        useUnk=TRUE
      }
    }
    #Are there marked no ID guys?
    useMarkednoID=FALSE
    if("G.marked.noID"%in%names(data)){
      if(!is.na(data$G.marked.noID[1])){
        G.marked.noID=data$G.marked.noID
        if(!is.matrix(G.marked.noID)){
          G.marked.noID=matrix(G.marked.noID)
        }
        if(ncol(G.marked.noID)!=ncat){
          stop("G.marked.noID needs ncat number of columns")
        }
        y.sight.marked.noID=data$y.sight.marked.noID
        useMarkednoID=TRUE
      }
    }

    #data checks
    if(length(dim(y.sight.marked))!=3){
      stop("dim(y.sight.marked) must be 3. Reduced to 2 during initialization")
    }
    if(length(dim(y.sight.unmarked))!=3){
      stop("dim(y.sight.unmarked) must be 3. Reduced to 2 during initialization")
    }
    if(useUnk){
      if(length(dim(y.sight.unk))!=3){
        stop("dim(y.sight.unk) must be 3. Reduced to 2 during initialization")
      }
    }
    if(useMarkednoID){
      if(length(dim(y.sight.marked.noID))!=3){
        stop("dim(y.sight.marked.noID) must be 3. Reduced to 2 during initialization")
      }
    }

    # if(is.na(nswap)){
    #   nswap=n.samp.latent/2
    #   warning("nswap not specified, using n.samp.latent/2")
    # }
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
    if(!(("psi1"%in%names(inits))|("psi2"%in%names(inits)))){
      stop("Must supply psi1 and psi2 inits.")
    }
    psi1<- inits$psi1
    psi2<- inits$psi2
    lam0=inits$lam0
    sigma<- inits$sigma
    gamma=inits$gamma
    if(!is.list(gamma)){
      stop("inits$gamma must be a list")
    }
    
    if(useUnk&!useMarkednoID){
      G.use=rbind(G.unmarked,G.unk)
      status=c(rep(2,nrow(G.unmarked)),rep(0,nrow(G.unk)))
      G.use=cbind(G.use,status)
      G.marked=cbind(G.marked,rep(1,nrow(G.marked)))
      ncat=ncat+1
      y.sight.latent=abind(y.sight.unmarked,y.sight.unk,along=1)
    }else if(!useUnk&useMarkednoID){
      G.use=rbind(G.unmarked,G.marked.noID)
      status=c(rep(2,nrow(G.unmarked)),rep(1,nrow(G.marked.noID)))
      G.use=cbind(G.use,status)
      G.marked=cbind(G.marked,rep(1,nrow(G.marked)))
      ncat=ncat+1
      y.sight.latent=abind(y.sight.unmarked,y.sight.marked.noID,along=1)
    }else if(useUnk&useMarkednoID){
      G.use=rbind(G.unmarked,G.unk,G.marked.noID)
      status=c(rep(2,nrow(G.unmarked)),rep(0,nrow(G.unk)),rep(1,nrow(G.marked.noID)))
      G.use=cbind(G.use,status)
      G.marked=cbind(G.marked,rep(1,nrow(G.marked)))
      ncat=ncat+1
      nIDcovs=c(nIDcovs,2)
      y.sight.latent=abind(y.sight.unmarked,y.sight.unk,y.sight.marked.noID,along=1)
    }else{
      G.use=G.unmarked
      y.sight.latent=y.sight.unmarked
    }
    n.samp.latent=nrow(y.sight.latent)
    if(is.na(nswap)){
      nswap=n.samp.latent/2
      warning("nswap not specified, using n.samp.latent/2")
    }
    
    #make constraints for data initialization
      constraints=matrix(1,nrow=n.samp.latent,ncol=n.samp.latent)
      for(i in 1:n.samp.latent){
        for(j in 1:n.samp.latent){
          guys1=which(G.use[i,]!=0)
          guys2=which(G.use[j,]!=0)
          comp=guys1[which(guys1%in%guys2)]
          if(any(G.use[i,comp]!=G.use[j,comp])){
            constraints[i,j]=0
          }
        }
      }
      #If bernoulli data, add constraints that prevent y.true[i,j,k]>1
      binconstraints=FALSE
      if(obstype=="bernoulli"){
        idx=which(y.sight.latent>0,arr.ind=TRUE)
        for(i in 1:n.samp.latent){
          for(j in 1:n.samp.latent){
            if(i!=j){
              if(all(idx[i,2:3]==idx[j,2:3])){
                constraints[i,j]=0 #can't combine samples from same trap and occasion in binomial model
                constraints[j,i]=0
                binconstraints=TRUE
              }
            }
          }
        }
      }
    
    #Build y.sight.true
    y.sight.marked.true=array(0,dim=c(M1,J,K))
    y.sight.unmarked.true=array(0,dim=c(M2,J,K))
    y.sight.marked.true[marked.guys,,]=y.sight.marked
    ID=matrix(NA,nrow=n.samp.latent,ncol=2) #second column: 1 indicates assigned to marked, 2 to unmarked
    idx=1
    for(i in 1:n.samp.latent){
      if(useMarkednoID){
        if(status[i]==1)next
      }
      if(idx>M2){
        stop("Need to raise M2 to initialize y.true")
      }
      traps=which(rowSums(y.sight.latent[i,,])>0)
      y.sight.unmarked.true2D=apply(y.sight.unmarked.true,c(1,2),sum)
      if(length(traps)==1){
        cand=which(y.sight.unmarked.true2D[,traps]>0)#guys caught at same traps
      }else{
        cand=which(rowSums(y.sight.unmarked.true2D[,traps])>0)#guys caught at same traps
      }
      if(length(cand)>0){
        if(length(cand)>1){#if more than 1 ID to match to, choose first one
          cand=cand[1]
        }
        #Check constraint matrix
        cands=which(ID[,1]%in%cand)#everyone assigned this ID
        if(all(constraints[i,cands]==1)){#focal consistent with all partials already assigned
          y.sight.unmarked.true[cand,,]=y.sight.unmarked.true[cand,,]+y.sight.latent[i,,]
          ID[i,]=c(cand,2)
        }else{#focal not consistent
          y.sight.unmarked.true[idx,,]=y.sight.latent[i,,]
          ID[i,]=c(idx,2)
          idx=idx+1
        }
      }else{#no assigned samples at this trap
        y.sight.unmarked.true[idx,,]=y.sight.latent[i,,]
        ID[i,]=c(idx,2)
        idx=idx+1
      }
    }
    if(useMarkednoID){#Need to initialize these guys to marked guys
      fix=which(status==1)
      meanloc=matrix(NA,nrow=n.marked,ncol=2)
      for(i in marked.guys){
        trap=which(rowSums(y.sight.marked[i,,])>0)
        locs2=matrix(0,nrow=0,ncol=2)
        if(length(trap)>0){
          locs2=rbind(locs2,X[trap,])
        }
        if(nrow(locs2)>1){
          meanloc[i,]=colMeans(locs2)
        }else if(nrow(locs2)>0){
          meanloc[i,]=locs2
        }
      }
      for(i in fix){
        trap=which(rowSums(y.sight.latent[i,,])>0)
        dists=sqrt((X[trap,1]-meanloc[,1])^2+(X[trap,2]-meanloc[,2])^2)
        ID[i,]=c(which(dists==min(dists,na.rm=TRUE))[1],1)
        y.sight.marked.true[ID[i],,]=y.sight.marked.true[ID[i],,]+y.sight.latent[i,,]
      }
    }
    
    #Check assignment consistency with constraints
    checkID=unique(ID)
    checkID=checkID[!checkID%in%marked.guys]
    for(i in 1:length(checkID)){
      idx=which(ID==checkID[i])
      if(!all(constraints[idx,idx]==1)){
        stop("ID initialized improperly")
      }
    }
   
    y.sight.marked.true=apply(y.sight.marked.true,c(1,2),sum)
    y.sight.unmarked.true=apply(y.sight.unmarked.true,c(1,2),sum)
    y.sight.latent=apply(y.sight.latent,c(1,2),sum)
    
    known.vector.marked=1*(rowSums(y.sight.marked.true>0))
    known.vector.unmarked=1*(rowSums(y.sight.unmarked.true>0))
    
    #Initialize zs
    z1=1*(known.vector.marked>0)
    z2=1*(known.vector.unmarked>0)
    add1=M1*(0.5-sum(z1)/M1)
    if(add1>0){
      z1[sample(which(z1==0),add1)]=1 #switch some uncaptured z's to 1.
    }
    add2=M2*(0.5-sum(z2)/M2)
    if(add2>0){
      z2[sample(which(z2==0),add2)]=1 #switch some uncaptured z's to 1.
    }
    
    #Optimize starting locations given where they are trapped.
    s1<- cbind(runif(M1,xlim[1],xlim[2]), runif(M1,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.sight.marked.true)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y.sight.marked.true[i,]>0,1:2],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s1[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s1[i,]<- trps
      }
    }
    s2<- cbind(runif(M2,xlim[1],xlim[2]), runif(M2,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.sight.unmarked.true)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y.sight.unmarked.true[i,]>0,1:2],ncol=2,byrow=FALSE)
      if(nrow(trps)>1){
        s2[i,]<- c(mean(trps[,1]),mean(trps[,2]))
      }else{
        s2[i,]<- trps
      }
    }
    
    #collapse unmarked data to 2D

    #Initialize G.true
    G.true.marked=matrix(0,nrow=M1,ncol=ncat)
    G.true.marked[marked.guys,]=G.marked
    G.true.unmarked=matrix(0,nrow=M2,ncol=ncat)
    for(i in unique(ID[,1])){
      idx=which(ID[,1]==i&ID[,2]==2)
      if(length(idx)==1){
        G.true.unmarked[i,]=G.use[idx,]
      }else{
        if(ncol(G.use)>1){
          G.true.unmarked[i,]=apply(G.use[idx,],2, max) #consensus
        }else{
          G.true.unmarked[i,]=max(G.use[idx,])
        }
      }
    }
    if(useUnk|useMarkednoID){#augmented guys are unmarked.
      if(max(ID)<M2){
        G.true.unmarked[,ncol(G.true.unmarked)]=2
      }
      unkguys=which(G.use[,ncol(G.use)]==0)
    }
    
    G.latent.marked=G.true.marked==0#Which genos can be updated?
    G.latent.unmarked=G.true.unmarked==0#Which genos can be updated?
    if(!(useUnk|useMarkednoID)){
      for(j in 1:(ncat)){
        fix=G.true.marked[,j]==0
        G.true.marked[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
        fix=G.true.unmarked[,j]==0
        G.true.unmarked[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
      }
    }else{
      for(j in 1:(ncat-1)){
        fix=G.true.marked[,j]==0
        G.true.marked[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
        fix=G.true.unmarked[,j]==0
        G.true.unmarked[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
      }
      #Split marked status back off
      Mark.obs=G.use[,ncat]
      # Mark.status=G.true[,ncat]
      ncat=ncat-1
      G.use=G.use[,1:ncat]
      G.true.marked=G.true.marked[,1:ncat]
      G.true.unmarked=G.true.unmarked[,1:ncat]
    }
    if(!is.matrix(G.use)){
      G.use=matrix(G.use,ncol=1)
    }
    if(!is.matrix(G.true.marked)){
      G.true.marked=matrix(G.true.marked,ncol=1)
    }
    if(!is.matrix(G.true.unmarked)){
      G.true.unmarked=matrix(G.true.unmarked,ncol=1)
    }
    # some objects to hold the MCMC output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=9)
    dimnames(out)<-list(NULL,c("lam0","sigma","n.m","n.um","N.m","N.um","N.all","psi1","psi2"))
    if(keepACs){
      s1xout<- s1yout<- z1out<-matrix(NA,nrow=nstore,ncol=M1)
      s2xout<- s2yout<- z2out<-matrix(NA,nrow=nstore,ncol=M2)
      IDout=array(NA,dim=c(nstore,nrow(ID),2))
    }
    iteridx=1 #for storing output not recorded every iteration
    if(keepGamma){
      gammaOut=vector("list",ncat)
      for(i in 1:ncat){
        gammaOut[[i]]=matrix(NA,nrow=nstore,ncol=nIDcovs[i])
        colnames(gammaOut[[i]])=paste("Lo",i,"G",1:nIDcovs[i],sep="")
      }
    }
    if(!is.na(data$locs[1])){
      uselocs=TRUE
      locs=data$locs
      telguys=which(rowSums(!is.na(locs[,,1]))>0)
      ll.tel=matrix(0,nrow=length(telguys),ncol=dim(locs)[2])
      for(i in telguys){
        ll.tel[i,]=dnorm(locs[i,,1],s1[i,1],sigma,log=TRUE)+dnorm(locs[i,,2],s1[i,2],sigma,log=TRUE)
      }
      ll.tel.cand=ll.tel
      #update starting locations using telemetry data
      for(i in telguys){
          s[i,]<- c(mean(locs[i,,1]),mean(locs[i,,2]))
      }
    }else{
      uselocs=FALSE
      telguys=c()
    }
    
    D1=e2dist(s1, X)
    D2=e2dist(s2, X)
    lamd.sight.marked<- lam0*exp(-D1*D1/(2*sigma*sigma))
    lamd.sight.unmarked<- lam0*exp(-D2*D2/(2*sigma*sigma))
    ll.y.sight.marked=array(0,dim=c(M1,J))
    ll.y.sight.unmarked=array(0,dim=c(M2,J))
    if(obstype=="bernoulli"){
      pd.sight.marked=1-exp(-lamd.sight.marked)
      pd.sight.unmarked=1-exp(-lamd.sight.unmarked)
      pd.sight.marked.cand=pd.sight.marked
      pd.sight.unmarked.cand=pd.sight.unmarked
      ll.y.sight.marked=dbinom(y.sight.marked.true,K,pd.sight.marked*z1,log=TRUE)
      ll.y.sight.unmarked=dbinom(y.sight.unmarked.true,K,pd.sight.unmarked*z2,log=TRUE)
    }else if(obstype=="poisson"){
      ll.y.sight.marked=dpois(y.sight.marked.true,K*lamd.sight.marked*z1,log=TRUE)
      ll.y.sight.unmarked=dpois(y.sight.unmarked.true,K*lamd.sight.unmarked*z2,log=TRUE)
    }
    lamd.sight.marked.cand=lamd.sight.marked
    lamd.sight.unmarked.cand=lamd.sight.unmarked
    ll.y.sight.marked.cand=ll.y.sight.marked
    ll.y.sight.unmarked.cand=ll.y.sight.unmarked
    
    for(iter in 1:niter){
      #Update both observation models
      if(obstype=="bernoulli"){
        #Update lam0
        llysightsum=sum(ll.y.sight.marked)+sum(ll.y.sight.unmarked)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.sight.marked.cand<- lam0.cand*exp(-D1*D1/(2*sigma*sigma))
          lamd.sight.unmarked.cand<- lam0.cand*exp(-D2*D2/(2*sigma*sigma))
          pd.sight.marked.cand=1-exp(-lamd.sight.marked.cand)
          pd.sight.unmarked.cand=1-exp(-lamd.sight.unmarked.cand)
          ll.y.sight.marked.cand= dbinom(y.sight.marked.true,K,pd.sight.marked.cand*z1,log=TRUE)
          ll.y.sight.unmarked.cand= dbinom(y.sight.unmarked.true,K,pd.sight.unmarked.cand*z2,log=TRUE)
          llysightcandsum=sum(ll.y.sight.marked.cand)+sum(ll.y.sight.unmarked.cand)
          if(runif(1) < exp(llysightcandsum-llysightsum)){
            lam0<- lam0.cand
            lamd.sight.marked=lamd.sight.marked.cand
            lamd.sight.unmarked=lamd.sight.unmarked.cand
            pd.sight.marked=pd.sight.marked.cand
            pd.sight.unmarked=pd.sight.unmarked.cand
            ll.y.sight.marked=ll.y.sight.marked.cand
            ll.y.sight.unmarked=ll.y.sight.unmarked.cand
            llysightsum=llysightcandsum
          }
        }
        #Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.sight.marked.cand<- lam0*exp(-D1*D1/(2*sigma.cand*sigma.cand))
          lamd.sight.unmarked.cand<- lam0*exp(-D2*D2/(2*sigma.cand*sigma.cand))
          pd.sight.marked.cand=1-exp(-lamd.sight.marked.cand)
          pd.sight.unmarked.cand=1-exp(-lamd.sight.unmarked.cand)
          ll.y.sight.marked.cand= dbinom(y.sight.marked.true,K,pd.sight.marked.cand*z1,log=TRUE)
          ll.y.sight.unmarked.cand= dbinom(y.sight.unmarked.true,K,pd.sight.unmarked.cand*z2,log=TRUE)
          llysightcandsum=sum(ll.y.sight.marked.cand)+sum(ll.y.sight.unmarked.cand)
          if(uselocs){
            for(i in telguys){
              ll.tel.cand[i,]=dnorm(locs[i,,1],s[i,1],sigma.cand,log=TRUE)+dnorm(locs[i,,2],s[i,2],sigma.cand,log=TRUE)
            }
          }else{
            ll.tel.cand=ll.tel=0
          }
          if(runif(1) < exp((llysightcandsum+sum(ll.tel.cand))-(llysightsum+sum(ll.tel)))){
            sigma<- sigma.cand
            lamd.sight.marked=lamd.sight.marked.cand
            lamd.sight.unmarked=lamd.sight.unmarked.cand
            pd.sight.marked=pd.sight.marked.cand
            pd.sight.unmarked=pd.sight.unmarked.cand
            ll.y.sight.marked=ll.y.sight.marked.cand
            ll.y.sight.unmarked=ll.y.sight.unmarked.cand
            ll.tel=ll.tel.cand
          }
        }
      }else{#poisson
        #Update lam0
        llysightsum=sum(ll.y.sight.marked)+sum(ll.y.sight.unmarked)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.sight.marked.cand<- lam0.cand*exp(-D1*D1/(2*sigma*sigma))
          lamd.sight.unmarked.cand<- lam0.cand*exp(-D2*D2/(2*sigma*sigma))
          ll.y.sight.marked.cand= dpois(y.sight.marked.true,K*lamd.sight.marked.cand*z1,log=TRUE)
          ll.y.sight.unmarked.cand= dpois(y.sight.unmarked.true,K*lamd.sight.unmarked.cand*z2,log=TRUE)
          llysightcandsum=sum(ll.y.sight.marked.cand)+sum(ll.y.sight.unmarked.cand)
          if(runif(1) < exp(llysightcandsum-llysightsum)){
            lam0<- lam0.cand
            lamd.sight.marked=lamd.sight.marked.cand
            lamd.sight.unmarked=lamd.sight.unmarked.cand
            ll.y.sight.marked=ll.y.sight.marked.cand
            ll.y.sight.unmarked=ll.y.sight.unmarked.cand
            llysightsum=llysightcandsum
          }
        }
        # Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.sight.marked.cand<- lam0*exp(-D1*D1/(2*sigma.cand*sigma.cand))
          lamd.sight.unmarked.cand<- lam0*exp(-D2*D2/(2*sigma.cand*sigma.cand))
          ll.y.sight.marked.cand= dpois(y.sight.marked.true,K*lamd.sight.marked.cand*z1,log=TRUE)
          ll.y.sight.unmarked.cand= dpois(y.sight.unmarked.true,K*lamd.sight.unmarked.cand*z2,log=TRUE)
          llysightcandsum=sum(ll.y.sight.marked.cand)+sum(ll.y.sight.unmarked.cand)
          if(uselocs){
            for(i in telguys){
              ll.tel.cand[i,]=dnorm(locs[i,,1],s[i,1],sigma.cand,log=TRUE)+dnorm(locs[i,,2],s[i,2],sigma.cand,log=TRUE)
            }
          }else{
            ll.tel.cand=ll.tel=0
          }
          if(runif(1) < exp((llysightcandsum+sum(ll.tel.cand))-(llysightsum+sum(ll.tel)))){  
            sigma<- sigma.cand
            lamd.sight.marked=lamd.sight.marked.cand
            lamd.sight.unmarked=lamd.sight.unmarked.cand
            ll.y.sight.marked=ll.y.sight.marked.cand
            ll.y.sight.unmarked=ll.y.sight.unmarked.cand
            ll.tel=ll.tel.cand
          }
        }
      }
      
      # ID update
      if(IDup=="Gibbs"){
        #Update y.sight.true from full conditional canceling out inconsistent combos with constraints.
        up=sample(1:n.samp.latent,nswap,replace=FALSE)
        for(l in up){
          nj=which(y.sight.latent[l,]>0)
          #Can only swap if IDcovs match
          idx=which(G.use[l,]!=0)
          if(length(idx)>1){#multiple loci observed
            possible.M=which(z1==1&apply(G.true.marked[,idx],1,function(x){all(x==G.use[l,idx])}))
            possible.UM=which(z2==1&apply(G.true.unmarked[,idx],1,function(x){all(x==G.use[l,idx])}))
          }else if(length(idx)==1){#single loci observed
            possible.M=which(z1==1&G.true.marked[,idx]==G.use[l,idx])
            possible.UM=which(z2==1&G.true.unmarked[,idx]==G.use[l,idx])
          }else{#fully latent G.obs
            possible.M=which(z1==1)#Can match anyone
            possible.UM=which(z2==1)#Can match anyone
          }
          if(!(useUnk|useMarkednoID)){#mark status exclusions handled through G.true
            possible.M=c() #Can't swap to a marked guy
          }else{
            if(Mark.obs[l]==2){#This is an unmarked sample
              possible.M=c()#Can't swap to a marked guy
            }
            if(Mark.obs[l]==1){#This is a marked sample
              possible.UM=c()#Can't swap to an unmarked guy
            }
          }
          if(length(c(possible.M,possible.UM))==0)next
          njprobs.M=lamd.sight.marked[,nj]
          njprobs.M[setdiff(1:M1,possible.M)]=0
          njprobs.UM=lamd.sight.unmarked[,nj]
          njprobs.UM[setdiff(1:M2,possible.UM)]=0
          njprobs=c(njprobs.M,njprobs.UM)
          njprobs=njprobs/sum(njprobs)
          # if(length(possible.M)==0){
          #   samptype="unmarked"
          # }else if(length(possible.UM)){
          #   samptype="marked"
          # }else{
          #   samptype="unk"
          # }
          
          newID=sample(1:(M1+M2),1,prob=njprobs)
          if(newID<=M1){
            newID=c(newID,1)
          }else{
            newID=c(newID-M1,2)
          }
      
          if(!all(ID[l,]==newID)){
            swapped=c(ID[l,1],newID[1])
            if(ID[l,2]==1&newID[2]==1){
              samptype="M2M"
            }else if(ID[l,2]==1&newID[2]==2){
              samptype="M2UM"
            }else if(ID[l,2]==2&newID[2]==2){
              samptype="UM2UM"
            }else if(ID[l,2]==2&newID[2]==1){
              samptype="UM2M"
            }
            if(samptype=="M2M"){#marked guy
              y.sight.marked.true[ID[l,1],]=y.sight.marked.true[ID[l,1],]-y.sight.latent[l,]
              y.sight.marked.true[newID[1],]=y.sight.marked.true[newID[1],]+y.sight.latent[l,]
              if(obstype=="bernoulli"){
                ll.y.marked.sight[swapped,]= dbinom(y.sight.marked.true[swapped,],K,pd.sight.marked[swapped,],log=TRUE)
              }else{
                ll.y.sight.marked[swapped,]= dpois(y.sight.marked.true[swapped,],K*lamd.sight.marked[swapped,],log=TRUE)
              }
            }else if(samptype=="UM2UM"){#unmarked guy
              y.sight.unmarked.true[ID[l,1],]=y.sight.unmarked.true[ID[l,1],]-y.sight.latent[l,]
              y.sight.unmarked.true[newID[1],]=y.sight.unmarked.true[newID[1],]+y.sight.latent[l,]
              if(obstype=="bernoulli"){
                ll.y.unmarked.sight[swapped,]= dbinom(y.sight.unmarked.true[swapped,],K,pd.sight.unmarked[swapped,],log=TRUE)
              }else{
                ll.y.sight.unmarked[swapped,]= dpois(y.sight.unmarked.true[swapped,],K*lamd.sight.unmarked[swapped,],log=TRUE)
              }
            }else if(samptype=="M2UM"){
              y.sight.marked.true[ID[l,1],]=y.sight.marked.true[ID[l,1],]-y.sight.latent[l,]
              y.sight.unmarked.true[newID[1],]=y.sight.unmarked.true[newID[1],]+y.sight.latent[l,]
              if(obstype=="bernoulli"){
                ll.y.marked.sight[swapped[1],]= dbinom(y.sight.marked.true[swapped[1],],K,pd.sight.marked[swapped[1],],log=TRUE)
                ll.y.unmarked.sight[swapped[2],]= dbinom(y.sight.unmarked.true[swapped[2],],K,pd.sight.unmarked[swapped[2],],log=TRUE)
              }else{
                ll.y.sight.marked[swapped[1],]= dpois(y.sight.marked.true[swapped[1],],K*lamd.sight.marked[swapped[1],],log=TRUE)
                ll.y.sight.unmarked[swapped[2],]= dpois(y.sight.unmarked.true[swapped[2],],K*lamd.sight.unmarked[swapped[2],],log=TRUE)
              }
            }else if(samptype=="UM2M"){
              y.sight.unmarked.true[ID[l,1],]=y.sight.unmarked.true[ID[l,1],]-y.sight.latent[l,]
              y.sight.marked.true[newID[1],]=y.sight.marked.true[newID[1],]+y.sight.latent[l,]
              if(obstype=="bernoulli"){
                ll.y.marked.sight[swapped[2],]= dbinom(y.sight.marked.true[swapped[2],],K,pd.sight.marked[swapped[2],],log=TRUE)
                ll.y.unmarked.sight[swapped[1],]= dbinom(y.sight.unmarked.true[swapped[1],],K,pd.sight.unmarked[swapped[1],],log=TRUE)
              }else{
                ll.y.sight.marked[swapped[2],]= dpois(y.sight.marked.true[swapped[2],],K*lamd.sight.marked[swapped[2],],log=TRUE)
                ll.y.sight.unmarked[swapped[1],]= dpois(y.sight.unmarked.true[swapped[1],],K*lamd.sight.unmarked[swapped[1],],log=TRUE)
              }
            }
            ID[l,]=newID
          }
        }
      }else{
        
        stop("MH for unknown number of marks not yet coded.")
        #### Fix this for unknown # marks
# 
#         up=sample(1:n.samp.latent,nswap,replace=FALSE)
#         y.sight.cand=y.sight.true
#         for(l in up){
#           #find legal guys to swap with. z=1 and consistent constraints
#           nj=which(y.sight.latent[l,]>0)
#           #Can only swap if IDcovs match
#           idX=which(G.use[l,]!=0)
#           if(length(idX)>1){#multiple loci observed
#             possible=which(z==1&apply(G.true[,idX],1,function(x){all(x==G.use[l,idX])}))
#           }else if(length(idX)==1){#single loci observed
#             possible=which(z==1&G.true[,idX]==G.use[l,idX])
#           }else{#fully latent G.obs
#             possible=which(z==1)#Can match anyone
#           }
#           if(!(useUnk|useMarkednoID)){#mark status exclusions handled through G.true
#             possible=possible[possible>n.marked]#Can't swap to a marked guy
#           }else{
#             if(Mark.obs[l]==2){#This is an unmarked sample
#               possible=possible[possible>n.marked]#Can't swap to a marked guy
#             }
#             if(Mark.obs[l]==1){#This is a marked sample
#               possible=possible[possible<=n.marked]#Can't swap to an unmarked guy
#             }
#           }
#           if(binconstraints){#can't have a y[i,j,k]>1
#             legal=rep(TRUE,length(possible))
#             for(i in 1:length(possible)){
#               check=which(ID==possible[i])#Who else is currently assigned this possible new ID?
#               if(length(check)>0){#if false, no samples assigned to this guy and legal stays true
#                 if(any(constraints[l,check]==0)){#if any members of the possible cluster are inconsistent with sample, illegal move
#                   legal[i]=FALSE
#                 }
#               }
#             }
#             possible=possible[legal]
#           }
#           if(length(possible)==0)next
#           njprobs=lamd.sight[,nj]
#           njprobs[setdiff(1:M,possible)]=0
#           njprobs=njprobs/sum(njprobs)
#           newID=ID
#           newID[l]=sample(1:M,1,prob=njprobs)
#           if(ID[l]==newID[l])next
# 
#           swapped=c(ID[l],newID[l])#order swap.out then swap.in
#           propprob=njprobs[swapped[2]]
#           backprob=njprobs[swapped[1]]
#           focalprob=1/n.samp.latent
#           focalbackprob=1/length(possible)
#           #update y.true
#           y.sight.cand[ID[l],]=y.sight.true[ID[l],]-y.sight.latent[l,]
#           y.sight.cand[newID[l],]=y.sight.true[newID[l],]+y.sight.latent[l,]
#           ##update ll.y
#           if(obstype=="poisson"){
#             ll.y.sight.cand[swapped,]=dpois(y.sight.cand[swapped,],K*lamd.sight[swapped,],log=TRUE)
#           }else{
#             ll.y.sight.cand[swapped,]=dbinom(y.sight.cand[swapped,],K,pd.sight[swapped,],log=TRUE)
#           }
#           if(runif(1)<exp(sum(ll.y.sight.cand[swapped,])-sum(ll.y.sight[swapped,]))*
#              (backprob/propprob)*(focalbackprob/focalprob)){
#             y.sight.true[swapped,]=y.sight.cand[swapped,]
#             ll.y.sight[swapped,]=ll.y.sight.cand[swapped,]
#             ID[l]=newID[l]
#           }
#         }
      }
      #update known.vector
      known.vector.marked=1*(rowSums(y.sight.marked.true>0))
      known.vector.unmarked=1*(rowSums(y.sight.unmarked.true>0))
      
      #Update G.true.marked
      G.true.marked.tmp=matrix(0, nrow=M1,ncol=ncat)
      G.true.marked.tmp[marked.guys,]=1
      for(i in unique(ID[ID[,2]==1,1])){
        idX=which(ID[,1]==i&ID[,2]==1)
        if(length(idX)==1){
          G.true.marked.tmp[i,]=G.use[idX,]
        }else{
          if(ncol(G.use)>1){
            G.true.marked.tmp[i,]=apply(G.use[idX,],2, max) #consensus
          }else{
            G.true.marked.tmp[i,]=max(G.use[idX,]) #consensus
          }
        }
      }
      G.latent.marked=G.true.marked.tmp==0
      for(j in 1:ncat){
        swap=G.latent.marked[,j]
        G.true.marked[swap,j]=sample(IDcovs[[j]],sum(swap),replace=TRUE,prob=gamma[[j]])
      }
      
      #Update G.true.unmarked
      G.true.unmarked.tmp=matrix(0, nrow=M2,ncol=ncat)
      for(i in unique(ID[ID[,2]==2,1])){
        idX=which(ID[,1]==i&ID[,2]==2)
        if(length(idX)==1){
          G.true.unmarked.tmp[i,]=G.use[idX,]
        }else{
          if(ncol(G.use)>1){
            G.true.unmarked.tmp[i,]=apply(G.use[idX,],2, max) #consensus
          }else{
            G.true.unmarked.tmp[i,]=max(G.use[idX,]) #consensus
          }
        }
      }
      G.latent.unmarked=G.true.unmarked.tmp==0
      for(j in 1:ncat){
        swap=G.latent.unmarked[,j]
        G.true.unmarked[swap,j]=sample(IDcovs[[j]],sum(swap),replace=TRUE,prob=gamma[[j]])
      }
      
      #update genotype frequencies
      for(j in 1:ncat){
        x=rep(NA,nIDcovs[[j]])
        for(k in 1:nIDcovs[[j]]){
          x[k]=sum(G.true.marked[z1==1,j]==k)+sum(G.true.unmarked[z2==1,j]==k)#genotype freqs in pop
        }
        gam=rgamma(rep(1,nIDcovs[[j]]),1+x)
        gamma[[j]]=gam/sum(gam)
      }
      
      ## probability of not being captured in a trap AT ALL by either method
      if(obstype=="poisson"){
        pd.sight.marked=1-exp(-lamd.sight.marked)
        pd.sight.unmarked=1-exp(-lamd.sight.unmarked)
      }
      pbar.sight.marked=(1-pd.sight.marked)^K
      pbar.sight.unmarked=(1-pd.sight.unmarked)^K
      prob0.sight.marked= exp(rowSums(log(pbar.sight.marked)))
      prob0.sight.unmarked= exp(rowSums(log(pbar.sight.unmarked)))
      fc.marked<- prob0.sight.marked*psi1/(prob0.sight.marked*psi1 + 1-psi1)
      fc.unmarked<- prob0.sight.unmarked*psi2/(prob0.sight.unmarked*psi2 + 1-psi2)
      z1[known.vector.marked==0]<- rbinom(sum(known.vector.marked ==0), 1, fc.marked[known.vector.marked==0])
      z2[known.vector.unmarked==0]<- rbinom(sum(known.vector.unmarked ==0), 1, fc.unmarked[known.vector.unmarked==0])
      if(obstype=="bernoulli"){
        ll.y.sight.marked= dbinom(y.sight.marked.true,K,pd.sight.marked*z1,log=TRUE)
        ll.y.sight.unmarked= dbinom(y.sight.unmarked.true,K,pd.sight.unmarked*z2,log=TRUE)
      }else{
        ll.y.sight.marked= dpois(y.sight.marked.true,K*lamd.sight.marked*z1,log=TRUE)
        ll.y.sight.unmarked= dpois(y.sight.unmarked.true,K*lamd.sight.unmarked*z2,log=TRUE)
      }
      psi1=rbeta(1,1+sum(z1),1+M1-sum(z1))
      psi2=rbeta(1,1+sum(z2),1+M2-sum(z2))
      
      ## Now we have to update the activity centers
      for (i in 1:M1) {#marked guys first
        if(i%in%telguys){
          Scand <- c(rnorm(1, s1[i, 1], proppars$st), rnorm(1, s1[i, 2], proppars$st))
        }else{
          Scand <- c(rnorm(1, s1[i, 1], proppars$s), rnorm(1, s1[i, 2], proppars$s))
        }
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamd.sight.marked.cand[i,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          if(obstype=="bernoulli"){
            pd.sight.marked.cand[i,]=1-exp(-lamd.sight.marked.cand[i,])
            ll.y.sight.marked.cand[i,]= dbinom(y.sight.marked.true[i,],K,pd.sight.marked.cand[i,]*z1[i],log=TRUE)
            if(uselocs&(i%in%telguys)){
              ll.tel.cand[i,]=dnorm(locs[i,,1],Scand[1],sigma,log=TRUE)+dnorm(locs[i,,2],Scand[2],sigma,log=TRUE)
              if (runif(1) < exp((sum(ll.y.sight.marked.cand[i,])+sum(ll.tel.cand[i,])) -
                                 (sum(ll.y.sight.marked[i,])+sum(ll.tel[i,])))) {
                s1[i,]=Scand
                D1[i,]=dtmp
                lamd.sight.marked[i,]=lamd.sight.cand.marked[i,]
                pd.sight.marked[i,]=pd.sight.marked.cand[i,]
                ll.y.sight.marked[i,]=ll.y.sight.marked.cand[i,]            
                ll.tel[i,]=ll.tel.cand[i,]
              }
            }else{
              if (runif(1) < exp(sum(ll.y.sight.marked.cand[i,]) -sum(ll.y.sight.marked[i,]))) {
                s1[i,]=Scand
                D1[i,]=dtmp
                lamd.sight.marked[i,]=lamd.sight.marked.cand[i,]
                pd.sight.marked[i,]=pd.sight.marked.cand[i,]
                ll.y.sight.marked[i,]=ll.y.sight.marked.cand[i,]            
              }
            }
          }else{#poisson
            ll.y.sight.marked.cand[i,]= dpois(y.sight.marked.true[i,],K*lamd.sight.marked.cand[i,]*z1[i],log=TRUE)
            if(uselocs&(i%in%telguys)){
              ll.tel.cand[i,]=dnorm(locs[i,,1],Scand[1],sigma,log=TRUE)+dnorm(locs[i,,2],Scand[2],sigma,log=TRUE)
              if (runif(1) < exp((sum(ll.y.sight.marked.cand[i,])+sum(ll.tel.cand[i,])) -
                                 (sum(ll.y.sight.marked[i,])+sum(ll.tel[i,])))) {
                s1[i,]=Scand
                D1[i,]=dtmp
                lamd.sight.marked[i,]=lamd.sight.marked.cand[i,]
                ll.y.sight.marked[i,]=ll.y.sight.marked.cand[i,]            
                ll.tel[i,]=ll.tel.cand[i,]
              }
            }else{
              if (runif(1) < exp(sum(ll.y.sight.marked.cand[i,]) -sum(ll.y.sight.marked[i,]))) {
                s1[i,]=Scand
                D1[i,]=dtmp
                lamd.sight.marked[i,]=lamd.sight.marked.cand[i,]
                ll.y.sight.marked[i,]=ll.y.sight.marked.cand[i,]            
              }
            }
          }
        }
      }
      ## Now unmarked guys
      for (i in 1:M2) {
        Scand <- c(rnorm(1, s2[i, 1], proppars$s), rnorm(1, s2[i, 2], proppars$s))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamd.sight.unmarked.cand[i,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          if(obstype=="bernoulli"){
            pd.sight.unmarked.cand[i,]=1-exp(-lamd.sight.unmarked.cand[i,])
            ll.y.sight.unmarked.cand[i,]= dbinom(y.sight.unmarked.true[i,],K,pd.sight.unmarked.cand[i,]*z2[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.sight.unmarked.cand[i,]) -sum(ll.y.sight.unmarked[i,]))) {
              s2[i,]=Scand
              D2[i,]=dtmp
              lamd.sight.unmarked[i,]=lamd.sight.unmarked.cand[i,]
              pd.sight.unmarked[i,]=pd.sight.unmarked.cand[i,]
              ll.y.sight.unmarked[i,]=ll.y.sight.unmarked.cand[i,]            
            }
          }else{#poisson
            ll.y.sight.unmarked.cand[i,]= dpois(y.sight.unmarked.true[i,],K*lamd.sight.unmarked.cand[i,]*z2[i],log=TRUE)
            
            if (runif(1) < exp(sum(ll.y.sight.unmarked.cand[i,]) -sum(ll.y.sight.unmarked[i,]))) {
              s2[i,]=Scand
              D2[i,]=dtmp
              lamd.sight.unmarked[i,]=lamd.sight.unmarked.cand[i,]
              ll.y.sight.unmarked[i,]=ll.y.sight.unmarked.cand[i,]            
            }
          }
          
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(keepACs){
          s1xout[iteridx,]<- s1[,1]
          s1yout[iteridx,]<- s1[,2]
          z1out[iteridx,]<- z1
          s2xout[iteridx,]<- s2[,1]
          s2yout[iteridx,]<- s2[,2]
          z2out[iteridx,]<- z2
          IDout[iteridx,,]=ID
        }
        if(keepGamma){
          for(k in 1:ncat){
            gammaOut[[k]][iteridx,]=gamma[[k]]
          }
        }
        NM=sum(z1)
        NUM=sum(z2)
        nM=length(unique(ID[ID[,2]==1,1])  )
        nUM=length(unique(ID[ID[,2]==2,1]))
        out[iteridx,]<- c(lam0,sigma,nM,nUM,NM,NUM,NM+NUM,psi1,psi2)
        iteridx=iteridx+1
      }
    }  # end of MCMC algorithm
    
    if(keepACs&keepGamma){
      list(out=out, s1xout=s1xout, s1yout=s1yout, z1out=z1out,s2xout=s2xout, s2yout=s2yout, z2out=z2out,IDout=IDout,gammaOut=gammaOut)
    }else if(keepACs&!keepGamma){
      list(out=out, s1xout=s1xout, s1yout=s1yout, z1out=z1out,s2xout=s2xout, s2yout=s2yout, z2out=z2out,IDout=IDout)
    }else if(!keepACs&keepGamma){
      list(out=out,gammaOut=gammaOut)
    }else{ 
      list(out=out)
    }
  }

