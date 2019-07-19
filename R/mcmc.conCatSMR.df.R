#' Fit the categorical spatial mark resight model with detection functions that vary by category level. Known
#' number of marked individuals.
#' @param data a data list as formatted by sim.conCatSMR.df(). See description for more details.
#' @param niter the number of MCMC iterations to perform
#' @param nburn the number of MCMC iterations to discard as burnin
#' @param nthin the MCMC thinning interval. Keep every nthin iterations.
#' @param M the level of data augmentation
#' @param inits a list of initial values for lam0, sigma, gamma, and psi. The list element for 
#' gamma is itself a list with ncat elements. The list elements for lam0
#' and sigma are also lists with starting values for each category level of the first categorical
#' identity covariate. See the example below.
#' @param obstype a character string indicating the observation model, "bernoulli" or "poisson".
#' @param nswap an integer indicating how many samples for which the latent identities
#' are updated on each iteration.
#' @param propars a list of proposal distribution tuning parameters for lam0, sigma, s, and st, for the
#' the activity centers of untelemetered and telemetered individuals, respectively. The tuning parameter
#' should be smaller for individuals with telemetry and increasingly so as the number of locations per
#' individual increases. The list elements for lam0 and sigma should also be lists with the same number
#' of levels specified in "inits".
#' @param storeLatent a logical indicator for whether or not the posteriors of the latent individual identities, z, and s are
#' stored and returned
#' @param storeGamma a logical indicator for whether or not the posteriors for gamma are stored and returned
#' @param IDup a character string indicating whether the latent identity update is done by Gibbs or Metropolis-
#' Hastings, "Gibbs", or "MH". For obstype="bernoulli", only "MH" is available because the full conditional is not known.
#' @description This function fits the conventional categorical spatial mark resight model when the number of
#' marked individuals is known and when the detection function parameters vary by one categorical
#' identity covariate. The distribution of marked individuals across the landscape is assumed to be
#' spatially uniform. This can be relaxed by modeling the marking process in the generalized spatial mark resight
#' samplers. An more simple version of this sampler that does not allow detection function parameters to vary by the
#' levels of one categorical covariate is located in mcmc.conCatSMR().
#' 
#' the data list should be formatted to match the list outputted by sim.conCatSMR(), but not all elements
#' of that object are necessary. y.sight.marked, y.sight.unmarked, G.marked, and G.unmarked are necessary
#' list elements. y.sight.x and G.x for x=unk and marke.noID are necessary if there are samples
#' of unknown marked status or samples from marked samples without individual identities.
#' 
#' An element "X", a matrix of trap coordinates, an element "K", the integer number of occasions, and
#' an element n.marked, the integer number of marked individuals are necessary.
#'
#' IDlist is a list containing elements ncat and IDcovs. ncat is an integer for the number
#' of categorical identity covariates and IDcovs is a list of length ncat with elements containing the
#' values each categorical identity covariate may take.
#' 
#' An element "locs", an n.marked x nloc x  2 array of telemetry locations is optional. This array can
#' have missing values if not all individuals have the same number of locations and the entry for individuals
#' with no telemetry should all be missing values (coded NA).
#' 
#' I will write a function to build the data object with "secr-like" input in the near future.
#' @author Ben Augustine
#' @examples
#' \dontrun{library(coda)
#' N=100
#' n.marked=40
#' lam0=c(0.1,0.3) #Enter same number of detection function parameters as nlevels[i]
#' sigma=c(0.7,0.5) #same number of lam0 and sigma parameters
#' K=20 #occasions
#' buff=3 #state space buffer
#' X<- expand.grid(3:11,3:11) #trapping array
#' pMarkID=c(1,1) #probability of observing marked status of marked and unmarked individuals
#' pID=1 #Probability marked individuals are identified
#' ncat=3  #number of ID categories
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each category level. Assume equal here.
#' nlevels=rep(2,ncat) #number of levels per IDcovs
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#category observation probabilities
#' tlocs=0 #telemetry locs/marked individual
#' obstype="poisson" #observation model, count or presence/absence?
#' marktype="premarked" #premarked or natural ID (marked individuals must be captured)?
#' data=sim.conCatSMR.df(N=N,n.marked=n.marked,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,ncat=ncat,
#'                          pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs,pID=pID,marktype=marktype)
#' 
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,psi=0.7)#start at simulated values
#' proppars=list(lam0=c(0.04,0.08),sigma=c(0.1,0.07),s=0.5,st=0.08)#proppar for each lam0 and sigma
#' 
#' M=175
#' storeLatent=TRUE
#' storeGamma=FALSE
#' niter=1500
#' nburn=0
#' nthin=1
#' IDup="Gibbs"
#' out=mcmc.conCatSMR.df(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,obstype=obstype,
#'                          proppars=proppars,storeLatent=TRUE,storeGamma=TRUE,IDup=IDup)
#' 
#' plot(mcmc(out$out))
#' 
#' 
#' 1-rejectionRate(mcmc(out$out)) #shoot for 0.2 - 0.4 for lam0 and sigma. If too low, raise proppar. If too high, lower proppar.
#' 1-rejectionRate(mcmc(out$sxout)) #activity center acceptance in x dimension. Shoot for min of 0.2
#' sum(rowSums(data$y.sight[(n.marked+1):N,,])>0) #true number of n.um
#' 
#' #Example of regular conventional SMR (no identity covariates)
#' #' N=100
#' n.marked=40
#' lam0=c(0.1,0.3) #Enter same number of detection function parameters as nlevels[i]
#' sigma=c(0.7,0.5) #same number of lam0 and sigma parameters
#' K=20 #occasions
#' buff=3 #state space buffer
#' X<- expand.grid(3:11,3:11) #trapping array
#' pMarkID=c(1,1) #probability of observing marked status of marked and unmarked individuals
#' pID=1 #Probability marked individuals are identified
#' ncat=1  #number of ID categories
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each category level. Assume equal here.
#' nlevels=rep(1,ncat) #number of levels per IDcovs
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#category observation probabilities
#' tlocs=0 #telemetry locs/marked individual
#' obstype="poisson" #observation model, count or presence/absence?
#' marktype="premarked" #premarked or natural ID (marked individuals must be captured)?
#' data=sim.conCatSMR.df(N=N,n.marked=n.marked,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype,ncat=ncat,
#'                          pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs,pID=pID,marktype=marktype)
#' 
#' inits=list(lam0=lam0,sigma=sigma,gamma=gamma,psi=0.7)#start at simulated values
#' proppars=list(lam0=c(0.04,0.08),sigma=c(0.1,0.07),s=0.5,st=0.08)#proppar for each lam0 and sigma
#' 
#' M=175
#' storeLatent=TRUE
#' storeGamma=FALSE
#' niter=1500
#' nburn=0
#' nthin=1
#' IDup="Gibbs"
#' out=mcmc.conCatSMR.df(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,obstype=obstype,
#'                          proppars=proppars,storeLatent=TRUE,storeGamma=TRUE,IDup=IDup)
#' 
#' plot(mcmc(out$out))
#' }
#' @export

mcmc.conCatSMR.df <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=NA,obstype="poisson",nswap=NA,
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           storeLatent=TRUE,storeGamma=TRUE,IDup="Gibbs",tf=NA){
    ###
    library(abind)
    y.sight.marked=data$y.sight.marked
    y.sight.unmarked=data$y.sight.unmarked
    X<-as.matrix(data$X)
    J<-nrow(X)
    K<- dim(y.sight.marked)[3]
    ncat=data$IDlist$ncat
    nallele=data$IDlist$nallele
    IDcovs=data$IDlist$IDcovs
    buff<- data$buff
    G.marked=data$G.marked
    n.marked=nrow(G.marked)
    if("G.unmarked"%in%names(data)){
      G.unmarked=data$G.unmarked
      if(length(dim(y.sight.unmarked))!=3){
        stop("dim(y.sight.unmarked) must be 3. Reduced to 2 during initialization")
      }
      useUM=TRUE
    }else{
      G.unmarked=matrix(0,nrow=0,ncol=ncat)
      useUM=FALSE
    }  
    if(!is.matrix(G.marked)){
      G.marked=matrix(G.marked)
    }
    if(!is.matrix(G.unmarked)){
      G.marked=matrix(G.unmarked)
    }
    if(!is.list(IDcovs)){
      stop("IDcovs must be a list")
    }
    nlevels=unlist(lapply(IDcovs,length))
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

    if(!IDup%in%c("MH","Gibbs")){
      stop("IDup must be MH or Gibbs")
    }
    if(IDup=="MH"){
      # stop("MH not implemented, yet")
    }
    if(obstype=="bernoulli"&IDup=="Gibbs"){
      stop("Must use MH IDup for bernoulli data")
    }
    #multiple df check
    ndf=length(IDcovs[[1]])
    if(length(inits$lam0)!=ndf|length(inits$sigma)!=ndf){
      stop("Must provide same number of detection function inits as category levels for first covariate")
    }
    if(length(proppars$lam0)!=ndf|length(proppars$sigma)!=ndf){
      stop("Must provide same number of detection function proposal parameters as category levels for first covariate")
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
    if(!any(is.na(tf))){
      if(any(tf>K)){
        stop("Some entries in tf are greater than K.")
      }
      if(is.null(dim(tf))){
        if(length(tf)!=J){
          stop("tf vector must be of length J.")
        }
        K2D=matrix(rep(tf,M),nrow=M,ncol=J,byrow=TRUE)
        warning("Since 1D tf entered, assuming all individuals exposed to equal capture")
      }else{
        if(!all(dim(tf)==c(M,J))){
          stop("tf must be dim M by J if tf varies by individual")
        }
        K2D=tf
        warning("Since 2D tf entered, assuming individual exposure to traps differ")
      }
    }else{
      tf=rep(K,J)
      K2D=matrix(rep(tf,M),nrow=M,ncol=J,byrow=TRUE)
    }
    
    
    ##pull out initial values
    psi<- inits$psi
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
      nlevels=c(nlevels,2)
      y.sight.latent=abind(y.sight.unmarked,y.sight.unk,y.sight.marked.noID,along=1)
    }else{
      G.use=G.unmarked
      y.sight.latent=y.sight.unmarked
    }
    n.samp.latent=nrow(y.sight.latent)
    if(is.na(nswap)){
      nswap=round(n.samp.latent/2)
      warning("nswap not specified, using round(n.samp.latent/2)")
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
        idx=t(apply(y.sight.latent,1,function(x){which(x>0,arr.ind=TRUE)}))
        for(i in 1:n.samp.latent){
          for(j in 1:n.samp.latent){
            if(i!=j){
              if(all(idx[i,1:2]==idx[j,1:2])){
                constraints[i,j]=0 #can't combine samples from same trap and occasion in binomial model
                constraints[j,i]=0
                binconstraints=TRUE
              }
            }
          }
        }
      }
    
    
    #Build y.sight.true
    y.sight.true=array(0,dim=c(M,J,K))
    y.sight.true[1:n.marked,,]=y.sight.marked
    ID=rep(NA,n.samp.latent)
    idx=n.marked+1
    for(i in 1:n.samp.latent){
      if(useMarkednoID){
        if(status[i]==1)next
      }
      if(idx>M){
        stop("Need to raise M to initialize y.true")
      }
      traps=which(rowSums(y.sight.latent[i,,])>0)
      y.sight.true2D=apply(y.sight.true,c(1,2),sum)
      if(length(traps)==1){
        cand=which(y.sight.true2D[,traps]>0)#guys caught at same traps
      }else{
        cand=which(rowSums(y.sight.true2D[,traps])>0)#guys caught at same traps
      }
      cand=cand[cand>n.marked]
      if(length(cand)>0){
        if(length(cand)>1){#if more than 1 ID to match to, choose first one
          cand=cand[1]
        }
        #Check constraint matrix
        cands=which(ID%in%cand)#everyone assigned this ID
        if(all(constraints[i,cands]==1)){#focal consistent with all partials already assigned
          y.sight.true[cand,,]=y.sight.true[cand,,]+y.sight.latent[i,,]
          ID[i]=cand
        }else{#focal not consistent
          y.sight.true[idx,,]=y.sight.latent[i,,]
          ID[i]=idx
          idx=idx+1
        }
      }else{#no assigned samples at this trap
        y.sight.true[idx,,]=y.sight.latent[i,,]
        ID[i]=idx
        idx=idx+1
      }
    }
    if(useMarkednoID){#Need to initialize these guys to marked guys
      fix=which(status==1)
      meanloc=matrix(NA,nrow=n.marked,ncol=2)
      for(i in 1:n.marked){
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
        ID[i]=which(dists==min(dists,na.rm=TRUE))[1]
        y.sight.true[ID[i],,]=y.sight.true[ID[i],,]+y.sight.latent[i,,]
      }
    }
    if(binconstraints){
      if(any(y.sight.true>1))stop("bernoulli data not initialized correctly")
    }
    
    #Check assignment consistency with constraints
    checkID=unique(ID)
    checkID=checkID[checkID>n.marked]
    for(i in 1:length(checkID)){
      idx=which(ID==checkID[i])
      if(!all(constraints[idx,idx]==1)){
        stop("ID initialized improperly")
      }
    }
   
    y.sight.true=apply(y.sight.true,c(1,2),sum)
    known.vector=c(rep(1,n.marked),rep(0,M-n.marked))
    known.vector[(n.marked+1):M]=1*(rowSums(y.sight.true[(n.marked+1):M,])>0)
    
    #Initialize z
    z=1*(known.vector>0)
    add=M*(0.5-sum(z)/M)
    if(add>0){
      z[sample(which(z==0),add)]=1 #switch some uncaptured z's to 1.
    }
    unmarked=c(rep(FALSE,n.marked),rep(TRUE,M-n.marked))
    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.sight.true)>0) #switch for those actually caught
    for(i in idx){
      trps<- matrix(X[y.sight.true[i,]>0,1:2],ncol=2,byrow=FALSE)
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
    #collapse unmarked data to 2D
    y.sight.latent=apply(y.sight.latent,c(1,2),sum)
    
    #Initialize G.true
    G.true=matrix(0,nrow=M,ncol=ncat)
    G.true[1:n.marked,]=G.marked
    for(i in unique(ID)){
      idx=which(ID==i)
      if(length(idx)==1){
        G.true[i,]=G.use[idx,]
      }else{
        if(ncol(G.use)>1){
          G.true[i,]=apply(G.use[idx,],2, max) #consensus
        }else{
          G.true[i,]=max(G.use[idx,])
        }
      }
    }
    if(useUnk|useMarkednoID){#augmented guys are unmarked.
      if(max(ID)<M){
        G.true[(max(ID)+1):M,ncol(G.true)]=2
      }
      unkguys=which(G.use[,ncol(G.use)]==0)
    }
    
    G.latent=G.true==0#Which genos can be updated?
    if(!(useUnk|useMarkednoID)){
      for(j in 1:(ncat)){
        fix=G.true[,j]==0
        G.true[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
      }
    }else{
      for(j in 1:(ncat-1)){
        fix=G.true[,j]==0
        G.true[fix,j]=sample(IDcovs[[j]],sum(fix),replace=TRUE,prob=gamma[[j]])
      }
      #Split marked status back off
      Mark.obs=G.use[,ncat]
      # Mark.status=G.true[,ncat]
      ncat=ncat-1
      G.use=G.use[,1:ncat]
      G.true=G.true[,1:ncat]
    }
    if(!is.matrix(G.use)){
      G.use=matrix(G.use,ncol=1)
    }
    if(!is.matrix(G.true)){
      G.true=matrix(G.true,ncol=1)
    }
    # some objects to hold the MCMC output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=ndf*2+3)
    dimnames(out)<-list(NULL,c(paste("lam0",1:ndf),paste("sigma",1:ndf),"N","n.um","psi"))
    
    if(storeLatent){
      sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
      IDout=matrix(NA,nrow=nstore,ncol=length(ID))
    }
    iteridx=1 #for storing output not recorded every iteration
    if(storeGamma){
      gammaOut=vector("list",ncat)
      for(i in 1:ncat){
        gammaOut[[i]]=matrix(NA,nrow=nstore,ncol=nlevels[i])
        colnames(gammaOut[[i]])=paste("Lo",i,"G",1:nlevels[i],sep="")
      }
    }
    if(!is.na(data$locs[1])){
      uselocs=TRUE
      locs=data$locs
      telguys=which(rowSums(!is.na(locs[,,1]))>0)
      #update starting locations using telemetry data
      for(i in telguys){
          s[i,]<- c(mean(locs[i,,1],na.rm=TRUE),mean(locs[i,,2],na.rm=TRUE))
      }
      ll.tel=matrix(0,nrow=max(telguys),ncol=dim(locs)[2])
      for(i in telguys){
        ll.tel[i,]=dnorm(locs[i,,1],s[i,1],sigma[G.true[i,1]],log=TRUE)+dnorm(locs[i,,2],s[i,2],sigma[G.true[i,1]],log=TRUE)
      }
      ll.tel.cand=ll.tel
    }else{
      uselocs=FALSE
      telguys=c()
    }
    
    D=e2dist(s, X)
    lamd.sight<- lam0[1]*exp(-D*D/(2*sigma[1]*sigma[1]))
    for(i in 2:ndf){
      lamd.sight[G.true[,1]==i,]=lam0[i]*exp(-D[G.true[,1]==i,]^2/(2*sigma[i]*sigma[i]))
    }
    ll.y.sight=array(0,dim=c(M,J))
    if(obstype=="bernoulli"){
      pd.sight=1-exp(-lamd.sight)
      pd.sight.cand=pd.sight
      ll.y.sight=dbinom(y.sight.true,K2D,pd.sight*z,log=TRUE)
    }else if(obstype=="poisson"){
      ll.y.sight=dpois(y.sight.true,K2D*lamd.sight*z,log=TRUE)
    }
    lamd.sight.cand=lamd.sight
    ll.y.sight.cand=ll.y.sight
    if(!is.finite(sum(ll.y.sight)))stop("Starting likelihood not finite. 
                                  Try raising lam0 and/or sigma inits.")
    ll.cat=rep(0,M)
    for(i in 1:M){
      mn=rep(0,ndf)
      mn[G.true[i,1]]=1
      ll.cat[i]=dmultinom(mn,prob=gamma[[1]],log=TRUE)
    }
    ll.cat.cand=ll.cat
    
    for(iter in 1:niter){
      #Update both observation models
      llysightsum=sum(ll.y.sight)
      lamd.sight.cand=lamd.sight
      if(obstype=="bernoulli"){
        for(i in 1:ndf){
          #Update lam0
          lam0.cand<- rnorm(1,lam0[i],proppars$lam0[i])
          if(lam0.cand > 0){
            lamd.sight.cand[G.true[,1]==i,]= lam0.cand*exp(-D[G.true[,1]==i,]^2/(2*sigma[i]*sigma[i]))
            pd.sight.cand=1-exp(-lamd.sight.cand)
            ll.y.sight.cand= dbinom(y.sight.true,K2D,pd.sight.cand*z,log=TRUE)
            llysightcandsum=sum(ll.y.sight.cand)
            if(runif(1) < exp(llysightcandsum-llysightsum)){
              lam0[i]<- lam0.cand
              lamd.sight=lamd.sight.cand
              pd.sight=pd.sight.cand
              ll.y.sight=ll.y.sight.cand
              llysightsum=llysightcandsum
            }else{
              lamd.sight.cand=lamd.sight
            }
          }
          #Update sigma
          sigma.cand<- rnorm(1,sigma[i],proppars$sigma[i])
          if(sigma.cand > 0){
            lamd.sight.cand[G.true[,1]==i,]<- lam0[i]*exp(-D[G.true[,1]==i,]^2/(2*sigma.cand*sigma.cand))
            pd.sight.cand=1-exp(-lamd.sight.cand)
            ll.y.sight.cand= dbinom(y.sight.true,K2D,pd.sight.cand*z,log=TRUE)
            llysightcandsum=sum(ll.y.sight.cand)
            if(uselocs){
              ll.tel.cand=ll.tel
              for(j in telguys){
                if(G.true[j,1]==i){
                  ll.tel.cand[j,]=dnorm(locs[j,,1],s[j,1],sigma.cand,log=TRUE)+dnorm(locs[j,,2],s[j,2],sigma.cand,log=TRUE)
                }
              }
            }else{
              ll.tel.cand=ll.tel=0
            }
            if(runif(1) < exp((llysightcandsum+sum(ll.tel.cand,na.rm=TRUE))-(llysightsum+sum(ll.tel,na.rm=TRUE)))){
              sigma[i]<- sigma.cand
              lamd.sight=lamd.sight.cand
              pd.sight=pd.sight.cand
              ll.y.sight=ll.y.sight.cand
              ll.tel=ll.tel.cand
              llysightsum=llysightcandsum
            }else{
              lamd.sight.cand=lamd.sight
            }
          }
        }
      }else{#poisson
        #Update lam0
        for(i in 1:ndf){
          lam0.cand<- rnorm(1,lam0[i],proppars$lam0[i])
          if(lam0.cand > 0){
            lamd.sight.cand[G.true[,1]==i,]= lam0.cand*exp(-D[G.true[,1]==i,]^2/(2*sigma[i]*sigma[i]))
            ll.y.sight.cand= dpois(y.sight.true,K2D*lamd.sight.cand*z,log=TRUE)
            llysightcandsum=sum(ll.y.sight.cand)
            if(runif(1) < exp(llysightcandsum-llysightsum)){
              lam0[i]<- lam0.cand
              lamd.sight=lamd.sight.cand
              ll.y.sight=ll.y.sight.cand
              llysightsum=llysightcandsum
            }else{
              lamd.sight.cand=lamd.sight
            }
          }
          # Update sigma
          sigma.cand<- rnorm(1,sigma[i],proppars$sigma[i])
          if(sigma.cand > 0){
            lamd.sight.cand[G.true[,1]==i,]<- lam0[i]*exp(-D[G.true[,1]==i,]^2/(2*sigma.cand*sigma.cand))
            ll.y.sight.cand= dpois(y.sight.true,K2D*lamd.sight.cand*z,log=TRUE)
            llysightcandsum=sum(ll.y.sight.cand)
            if(uselocs){
              ll.tel.cand=ll.tel
              for(j in telguys){
                if(G.true[j,1]==i){
                  ll.tel.cand[j,]=dnorm(locs[j,,1],s[j,1],sigma.cand,log=TRUE)+dnorm(locs[j,,2],s[j,2],sigma.cand,log=TRUE)
                }
              }
            }else{
              ll.tel.cand=ll.tel=0
            }
            if(runif(1) < exp((llysightcandsum+sum(ll.tel.cand,na.rm=TRUE))-(llysightsum+sum(ll.tel,na.rm=TRUE)))){  
              sigma[i]<- sigma.cand
              lamd.sight=lamd.sight.cand
              ll.y.sight=ll.y.sight.cand
              llysightsum=llysightcandsum
              ll.tel=ll.tel.cand
            }else{
              lamd.sight.cand=lamd.sight
            }
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
            possible=which(z==1&apply(G.true[,idx],1,function(x){all(x==G.use[l,idx])}))
          }else if(length(idx)==1){#single loci observed
            possible=which(z==1&G.true[,idx]==G.use[l,idx])
          }else{#fully latent G.obs
            possible=which(z==1)#Can match anyone
          }
          if(!(useUnk|useMarkednoID)){#mark status exclusions handled through G.true
            possible=possible[possible>n.marked]#Can't swap to a marked guy
          }else{
            if(Mark.obs[l]==2){#This is an unmarked sample
              possible=possible[possible>n.marked]#Can't swap to a marked guy
            }
            if(Mark.obs[l]==1){#This is a marked sample
              possible=possible[possible<=n.marked]#Can't swap to an unmarked guy
            }
          }
          if(length(possible)==0)next
          njprobs=lamd.sight[,nj]
          njprobs[setdiff(1:M,possible)]=0
          njprobs=njprobs/sum(njprobs)
          newID=sample(1:M,1,prob=njprobs)
          if(ID[l]!=newID){
            swapped=c(ID[l],newID)
            #update y.true
            y.sight.true[ID[l],]=y.sight.true[ID[l],]-y.sight.latent[l,]
            y.sight.true[newID,]=y.sight.true[newID,]+y.sight.latent[l,]
            ID[l]=newID
            if(obstype=="bernoulli"){
              ll.y.sight[swapped,]= dbinom(y.sight.true[swapped,],K2D[swapped,],pd.sight[swapped,],log=TRUE)
            }else{
              ll.y.sight[swapped,]= dpois(y.sight.true[swapped,],K2D[swapped,]*lamd.sight[swapped,],log=TRUE)
            }
          }
        }
      }else{
        up=sample(1:n.samp.latent,nswap,replace=FALSE)
        y.sight.cand=y.sight.true
        for(l in up){
          #find legal guys to swap with. z=1 and consistent constraints
          nj=which(y.sight.latent[l,]>0)
          #Can only swap if IDcovs match
          idX=which(G.use[l,]!=0)
          if(length(idX)>1){#multiple loci observed
            possible=which(z==1&apply(G.true[,idX],1,function(x){all(x==G.use[l,idX])}))
          }else if(length(idX)==1){#single loci observed
            possible=which(z==1&G.true[,idX]==G.use[l,idX])
          }else{#fully latent G.obs
            possible=which(z==1)#Can match anyone
          }
          if(!(useUnk|useMarkednoID)){#mark status exclusions handled through G.true
            possible=possible[possible>n.marked]#Can't swap to a marked guy
          }else{
            if(Mark.obs[l]==2){#This is an unmarked sample
              possible=possible[possible>n.marked]#Can't swap to a marked guy
            }
            if(Mark.obs[l]==1){#This is a marked sample
              possible=possible[possible<=n.marked]#Can't swap to an unmarked guy
            }
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
          if(length(possible)==0)next
          njprobs=lamd.sight[,nj]
          njprobs[setdiff(1:M,possible)]=0
          njprobs=njprobs/sum(njprobs)
          newID=ID
          newID[l]=sample(1:M,1,prob=njprobs)
          if(ID[l]==newID[l])next

          swapped=c(ID[l],newID[l])#order swap.out then swap.in
          propprob=njprobs[swapped[2]]
          backprob=njprobs[swapped[1]]
          # focalprob=1/n.samp.latent
          # focalbackprob=1/length(possible)
          #update y.true
          y.sight.cand[ID[l],]=y.sight.true[ID[l],]-y.sight.latent[l,]
          y.sight.cand[newID[l],]=y.sight.true[newID[l],]+y.sight.latent[l,]
          focalprob=(sum(ID==ID[l])/n.samp.latent)*(y.sight.true[ID[l],nj]/sum(y.sight.true[ID[l],]))
          focalbackprob=(sum(newID==newID[l])/n.samp.latent)*(y.sight.cand[newID[l],nj]/sum(y.sight.cand[newID[l],]))
          ##update ll.y
          if(obstype=="poisson"){
            ll.y.sight.cand[swapped,]=dpois(y.sight.cand[swapped,],K2D[swapped,]*lamd.sight[swapped,],log=TRUE)
          }else{
            ll.y.sight.cand[swapped,]=dbinom(y.sight.cand[swapped,],K2D[swapped,],pd.sight[swapped,],log=TRUE)
          }
          if(runif(1)<exp(sum(ll.y.sight.cand[swapped,])-sum(ll.y.sight[swapped,]))*
             (backprob/propprob)*(focalbackprob/focalprob)){
            y.sight.true[swapped,]=y.sight.cand[swapped,]
            ll.y.sight[swapped,]=ll.y.sight.cand[swapped,]
            ID[l]=newID[l]
          }
        }
      }
      # #update known.vector and G.latent
      known.vector[(n.marked+1):M]=1*(rowSums(y.sight.true[(n.marked+1):M,])>0)
      G.true.tmp=matrix(0, nrow=M,ncol=ncat)
      G.true.tmp[1:n.marked,]=1
      for(i in unique(ID[ID>n.marked])){
        idX=which(ID==i)
        if(length(idX)==1){
          G.true.tmp[i,]=G.use[idX,]
        }else{
          if(ncol(G.use)>1){
            G.true.tmp[i,]=apply(G.use[idX,],2, max) #consensus
          }else{
            G.true.tmp[i,]=max(G.use[idX,]) #consensus
          }
        }
      }
      G.latent=G.true.tmp==0
      #update G.true with no df attached
      if(ncat>1){
        for(j in 2:ncat){
          swap=G.latent[,j]
          G.true[swap,j]=sample(IDcovs[[j]],sum(swap),replace=TRUE,prob=gamma[[j]])
        }
      }else{
        swap=G.latent[,1]
        G.true[swap,1]=sample(IDcovs[[1]],sum(swap),replace=TRUE,prob=gamma[[1]])
      }
      #update df loci
      idx=which(G.latent[,1])
      if(length(idx)>nswap){
        updf=sample(idx,nswap)
      }else{
        updf=idx
      }
      for(i in updf){
        G.cand=rmultinom(1,1,gamma[[1]])
        G.cand=which(G.cand==1)
        if(G.cand==G.true[i,1])next
        #update ll.y
        lamd.sight.cand[i,]<- lam0[G.cand]*exp(-D[i,]^2/(2*sigma[G.cand]*sigma[G.cand]))
        if(obstype=="bernoulli"){
          pd.sight.cand[i,]=1-exp(-lamd.sight.cand[i,])
          ll.y.sight.cand[i,]= dbinom(y.sight.true[i,],K2D[i,],pd.sight.cand[i,]*z[i],log=TRUE)
        }else{
          ll.y.sight.cand[i,]= dpois(y.sight.true[i,],K2D[i,]*lamd.sight.cand[i,]*z[i],log=TRUE)
        }
        mn=rep(0,ndf)
        mn[G.cand]=1
        ll.cat.cand[i]=dmultinom(mn,prob=gamma[[1]],log=TRUE)
        prop.for=gamma[[1]][G.cand]
        prop.back=gamma[[1]][G.true[i,1]]
        
        if(runif(1) < exp((ll.cat.cand[i]+sum(ll.y.sight.cand[i,]))-
                          (ll.cat[i]+sum(ll.y.sight[i,])))*(prop.back/prop.for)){
          G.true[i,1]=G.cand
          ll.cat[i]=ll.cat.cand[i]
          lamd.sight[i,]=lamd.sight.cand[i,]
          if(obstype=="bernoulli"){
            pd.sight[i,]=pd.sight.cand[i,]
          }
          ll.y.sight[i,]=ll.y.sight.cand[i,]
        }
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
      #update ll.cat
      for(i in 1:M){
        mn=rep(0,ndf)
        mn[G.true[i,1]]=1
        ll.cat[i]=dmultinom(mn,prob=gamma[[1]],log=TRUE)
      }
      ## probability of not being captured in a trap AT ALL by either method
      if(obstype=="poisson"){
        pd.sight=1-exp(-lamd.sight)
      }
      pbar.sight=(1-pd.sight)^K2D
      prob0.sight<- exp(rowSums(log(pbar.sight)))
      fc<- prob0.sight*psi/(prob0.sight*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      if(obstype=="bernoulli"){
        ll.y.sight= dbinom(y.sight.true,K2D,pd.sight*z,log=TRUE)
      }else{
        ll.y.sight= dpois(y.sight.true,K2D*lamd.sight*z,log=TRUE)
      }
      psi=rbeta(1,1+sum(z[unmarked]),1+M-n.marked-sum(z[unmarked]))
      
      ## Now we have to update the activity centers
      for (i in 1:M) {
        if(i%in%telguys){
          Scand <- c(rnorm(1, s[i, 1], proppars$st), rnorm(1, s[i, 2], proppars$st))
        }else{
          Scand <- c(rnorm(1, s[i, 1], proppars$s), rnorm(1, s[i, 2], proppars$s))
        }
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          dfi=G.true[i,1]
          lamd.sight.cand[i,]<- lam0[dfi]*exp(-dtmp*dtmp/(2*sigma[dfi]*sigma[dfi]))
          if(obstype=="bernoulli"){
            pd.sight.cand[i,]=1-exp(-lamd.sight.cand[i,])
            ll.y.sight.cand[i,]= dbinom(y.sight.true[i,],K2D[i,],pd.sight.cand[i,]*z[i],log=TRUE)
            if(uselocs&(i%in%telguys)){
              ll.tel.cand[i,]=dnorm(locs[i,,1],Scand[1],sigma[dfi],log=TRUE)+dnorm(locs[i,,2],Scand[2],sigma[dfi],log=TRUE)
              if (runif(1) < exp((sum(ll.y.sight.cand[i,])+sum(ll.tel.cand[i,],na.rm=TRUE)) -
                                 (sum(ll.y.sight[i,])+sum(ll.tel[i,],na.rm=TRUE)))) {
                s[i,]=Scand
                D[i,]=dtmp
                lamd.sight[i,]=lamd.sight.cand[i,]
                pd.sight[i,]=pd.sight.cand[i,]
                ll.y.sight[i,]=ll.y.sight.cand[i,]            
                ll.tel[i,]=ll.tel.cand[i,]
              }
            }else{
              if (runif(1) < exp(sum(ll.y.sight.cand[i,]) -sum(ll.y.sight[i,]))) {
                s[i,]=Scand
                D[i,]=dtmp
                lamd.sight[i,]=lamd.sight.cand[i,]
                pd.sight[i,]=pd.sight.cand[i,]
                ll.y.sight[i,]=ll.y.sight.cand[i,]            
              }
            }
          }else{#poisson
            ll.y.sight.cand[i,]= dpois(y.sight.true[i,],K2D[i,]*lamd.sight.cand[i,]*z[i],log=TRUE)
            if(uselocs&(i%in%telguys)){
              ll.tel.cand[i,]=dnorm(locs[i,,1],Scand[1],sigma[dfi],log=TRUE)+dnorm(locs[i,,2],Scand[2],sigma[dfi],log=TRUE)
              if (runif(1) < exp((sum(ll.y.sight.cand[i,])+sum(ll.tel.cand[i,],na.rm=TRUE)) -
                                 (sum(ll.y.sight[i,])+sum(ll.tel[i,],na.rm=TRUE)))) {
                s[i,]=Scand
                D[i,]=dtmp
                lamd.sight[i,]=lamd.sight.cand[i,]
                ll.y.sight[i,]=ll.y.sight.cand[i,]            
                ll.tel[i,]=ll.tel.cand[i,]
              }
            }else{
              if (runif(1) < exp(sum(ll.y.sight.cand[i,]) -sum(ll.y.sight[i,]))) {
                s[i,]=Scand
                D[i,]=dtmp
                lamd.sight[i,]=lamd.sight.cand[i,]
                ll.y.sight[i,]=ll.y.sight.cand[i,]            
              }
            }
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        if(storeLatent){
          sxout[iteridx,]<- s[,1]
          syout[iteridx,]<- s[,2]
          zout[iteridx,]<- z
          IDout[iteridx,]=ID
        }
        if(storeGamma){
          for(k in 1:ncat){
            gammaOut[[k]][iteridx,]=gamma[[k]]
          }
        }
        if(useUnk&(!useUM)){
          n=length(unique(ID))
        }else if(useUnk|useMarkednoID){
          n=length(unique(ID[ID>n.marked]))
        }else{
          n=length(unique(ID))
        }
        
        out[iteridx,]<- c(lam0,sigma,sum(z),n,psi)
        iteridx=iteridx+1
      }
    }  # end of MCMC algorithm
    
    if(storeLatent&storeGamma){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout,gammaOut=gammaOut)
    }else if(storeLatent&!storeGamma){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else if(!storeLatent&storeGamma){
      list(out=out,gammaOut=gammaOut)
    }else{ 
      list(out=out)
    }
  }

