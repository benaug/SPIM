#' Fit the generalized categorical spatial mark resight model allowing for activity center relocation
#' between the marking and sighting process
#' @param data a data list as formatted by sim.genCatSMR(). See description for more details.
#' @param niter the number of MCMC iterations to perform
#' @param nburn the number of MCMC iterations to discard as burnin
#' @param nthin the MCMC thinning interval. Keep every nthin iterations.
#' @param M the level of data augmentation
#' @param inits a list of initial values for lam0.mark,lam0.sight, sigma_d, sigma_p, gamma, and psi. The list element for 
#' gamma is itself a list with ncat elements. See the example below.
#' @param obstype a vector of length two indicating the observation model, "bernoulli" or "poisson", for the 
#' marking and sighting process
#' @param nswap an integer indicating how many samples for which the latent identities
#' are updated on each iteration.
#' @param propars a list of proposal distribution tuning parameters for lam0.mark, lam0.sight, sigma_d, sigma_p, s1, s2, and s2t, for the
#' the activity centers of the marking process, the sighting process untelemetered and sighting process telemetered individuals, respectively. The tuning parameter
#' should be smaller for individuals with telemetry and increasingly so as the number of locations per
#' individual increases
#' @param storeLatent a logical indicator for whether or not the posteriors of the latent individual identities, z, and s are
#' stored and returned
#' @param storeGamma a logical indicator for whether or not the posteriors for gamma are stored and returned
#' @param IDup a character string indicating whether the latent identity update is done by Gibbs or Metropolis-
#' Hastings, "Gibbs", or "MH". For obstype="bernoulli", only "MH" is available because the full conditional is not known.
#' @param tf1 a trap operation vector or matrix for the marking process. If exposure to capture does
#' not vary by indiviudal, tf1 should be a vector of length J1 indicating how many of the K1 occasions
#' each marking location was operational. If exposure to capture varies by individual or by trap and
#' individual, tf1 should be a matrix of dimension M x J1 indicating how many of the K1 occasions individual
#' i was exposed to at trap j. This allows known additions or removals during the marking process
#'  to be accounted for. Exposure for the n.marked+1 ... M uncaptured individuals should be the
#'  same as the number of occasions each trap was operational. We can't account for unknown
#'  additions and removals.
#' @param tf2 a trap operation vector or matrix for the sighting process. If exposure to capture does
#' not vary by indiviudal, tf1 should be a vector of length J2 indicating how many of the K2 occasions
#' each sighting location was operational. If exposure to capture varies by individual or by trap and
#' individual, tf2 should be a matrix of dimension M x J2 indicating how many of the K2 occasions individual
#' i was exposed to at trap j. This allows known additions or removals between the marking and
#' sighting processes and during the sighting process to be accounted for. Exposure for 
#' the n.marked+1 ... M uncaptured individuals should be the
#'  same as the number of occasions each trap was operational. We can't account for unknown
#'  additions and removals.
#' @description This function fits the generalized categorical spatial mark resight model that
#' allows the activity centers to move between marking and sighting process following a bivariate
#' normal Markov transition kernel.
#' Modelling the marking process relaxes the assumption that the distribution of marked 
#' individuals across the landscape is
#' spatially uniform. Allowing the activity centers to move between marking and sighting accommodates
#' animal movement, while keeping the marked individuals' activity centers near the marking process
#' capture locations. This is an intermediate model between having fixed activity centers and 
#' random activity center relocation between processes, which is unrealistic. A model with no 
#' activity center relocation is found in mcmc.genCatSMR() and that model with detection function 
#' parameters that vary by the identity category levels is found in mcmc.genCatSMR.df(). Mobile
#' activity centers and variable detection functions could be combined into one model.
#' 
#' the data list should be formatted to match the list outputted by sim.genCatSMR.move(), but not all elements
#' of that object are necessary. y.mark, y.sight.marked, y.sight.unmarked, G.marked, and G.unmarked are necessary
#' list elements. y.sight.x and G.x for x=unk and marke.noID are necessary if there are samples
#' of unknown marked status or samples from marked samples without individual identities.
#' 
#' An element "X1", a matrix of marking coordinates, an element "X2", a matrix of sighting coordinates, ,
#' an element "K1", the integer number of marking occasions, and an element "K2", the integer number of sighting occasions
#'  are necessary.
#' 
#' IDlist is a list containing elements ncat and IDcovs. ncat is an integer for the number
#' of categorical identity covariates and IDcovs is a list of length ncat with elements containing the
#' values each categorical identity covariate may take.
#' 
#' An element "locs", an n.marked x nloc x  2 array of telemetry locations is optional. This array can
#' have missing values if not all individuals have the same number of locations and the entry for individuals
#' with no telemetry should all be missing values (coded NA). These telemetry locations are around the sighting process activity centers only. Marking process telemetry
#' could be added. In the models with no activity center relocation, it does not matter.
#'
#'This version does not allow for interspersion of the marking and sighting process because
#'it does not make sense for the activity centers to temporally switch back and forth between
#'the marking and sighting process activity centers.
#'   
#' I will write a function to build the data object with "secr-like" input in the near future.
#' 
#' @author Ben Augustine
#' @examples
#' \dontrun{
#' #Using categorical identity covariates
#' N=50
#' lam0.mark=0.05
#' lam0.sight=0.25
#' sigma_d=0.50
#' sigma_p=1
#' K1=10
#' K2=10
#' buff=2
#' X1<- expand.grid(3:11,3:11)
#' X2<- expand.grid(3:11+0.5,3:11+0.5)
#' pMarkID=c(0.8,0.8) #probability of observing marked status of marked and unmarked guys
#' pID=0.8
#' obstype=c("bernoulli","poisson")
#' ncat=3  #number of loci
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each genotype. Assume equal for now
#' nlevels=rep(2,ncat) #number of IDcovs per loci
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#loci amplification/observation probabilities
#' tlocs=25
#' data=sim.genCatSMR.move(N=N,lam0.mark=lam0.mark,lam0.sight=lam0.sight,sigma_d=sigma_d,sigma_p=sigma_p,K1=K1,
#'                         K2=K2,X1=X1,X2=X2,buff=buff,obstype=obstype,ncat=ncat,
#'                         pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs)
#' 
#' inits=list(lam0.mark=lam0.mark,lam0.sight=lam0.sight,sigma_d=sigma_d,sigma_p=sigma_p,gamma=gamma,psi=0.7)
#' proppars=list(lam0.mark=0.05,lam0.sight=0.1,sigma_d=0.02,sigma_p=0.2,s1=0.5,s2=0.25,s2t=0.1)#poisson-poisson
#' M=100
#' storeLatent=TRUE
#' storeGamma=FALSE
#' niter=500
#' nburn=0
#' nthin=1
#' IDup="MH"
#' out=mcmc.genCatSMR.move(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,obstype=obstype,
#'                         proppars=proppars,storeLatent=TRUE,storeGamma=TRUE,IDup=IDup)
#' 
#' plot(mcmc(out$out))
#' 1-rejectionRate(mcmc(out$out))
#' 1-rejectionRate(mcmc(out$s1xout))
#' 1-rejectionRate(mcmc(out$s2xout))
#' length(unique(data$IDum)) #true number of unmarked individuals captured
#' 
#' ####regular generalized SMR with no identity covariates
#' N=50
#' lam0.mark=0.05
#' lam0.sight=0.25
#' sigma_d=0.50
#' sigma_p=1
#' K1=10
#' K2=10
#' buff=2
#' X1<- expand.grid(3:11,3:11)
#' X2<- expand.grid(3:11+0.5,3:11+0.5)
#' pMarkID=c(0.8,0.8) #probability of observing marked status of marked and unmarked guys
#' pID=0.8
#' obstype=c("bernoulli","poisson")
#' ncat=1  #just 1 covariate
#' gamma=IDcovs=vector("list",ncat) #population frequencies of each genotype. Assume equal for now
#' nlevels=rep(1,ncat) #just 1 value. We have simplified to regular generalized SMR.
#' for(i in 1:ncat){
#'   gamma[[i]]=rep(1/nlevels[i],nlevels[i])
#'   IDcovs[[i]]=1:nlevels[i]
#' }
#' pIDcat=rep(1,ncat)#loci amplification/observation probabilities
#' tlocs=25
#' data=sim.genCatSMR.move(N=N,lam0.mark=lam0.mark,lam0.sight=lam0.sight,sigma_d=sigma_d,sigma_p=sigma_p,K1=K1,
#'                         K2=K2,X1=X1,X2=X2,buff=buff,obstype=obstype,ncat=ncat,
#'                         pIDcat=pIDcat,gamma=gamma,IDcovs=IDcovs,pMarkID=pMarkID,tlocs=tlocs)
#' 
#' inits=list(lam0.mark=lam0.mark,lam0.sight=lam0.sight,sigma_d=sigma_d,sigma_p=sigma_p,gamma=gamma,psi=0.7)
#' proppars=list(lam0.mark=0.05,lam0.sight=0.1,sigma_d=0.02,sigma_p=0.2,s1=0.5,s2=0.25,s2t=0.1)#poisson-poisson
#' M=100
#' storeLatent=TRUE
#' storeGamma=FALSE
#' niter=500
#' nburn=0
#' nthin=1
#' IDup="MH"
#' out=mcmc.genCatSMR.move(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,obstype=obstype,
#'                         proppars=proppars,storeLatent=TRUE,storeGamma=TRUE,IDup=IDup)
#' 
#' plot(mcmc(out$out))
#'}
#' @export

mcmc.genCatSMR.move <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=NA,obstype=c("bernoulli","poisson"),nswap=NA,
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),
           storeLatent=TRUE,storeGamma=TRUE,IDup="Gibbs",tf1=NA,tf2=NA){
    if(any(data$markedS==0|data$markedS==2)){#capture order constraints
      mcmc.genCatSMR.moveb(data,niter=niter,nburn=nburn,nthin=nthin,M=M,inits=inits,
                      obstype=obstype,nswap=nswap,proppars=proppars,
                      storeLatent=storeLatent,storeGamma=storeGamma,IDup=IDup,tf1=tf1,tf2=tf2)
    }else{#no capture order constraints
      mcmc.genCatSMR.movea(data,niter=niter,nburn=nburn,nthin=nthin,M=M,inits=inits,
                      obstype=obstype,nswap=nswap,proppars=proppars,
                      storeLatent=storeLatent,storeGamma=storeGamma,IDup=IDup,tf1=tf1,tf2=tf2)
    }
  }

