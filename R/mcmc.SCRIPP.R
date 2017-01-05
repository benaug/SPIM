#' Run MCMC algorithm for SCR model with inhomogenous density and 1 categorical density covariate. Discrete state space.
#' @param data a list produced by simSCRIPP or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam01=lam01,lam02=lam02,sigma=sigma)
#' @param proppars a list of tuning parameters for the proposal distributions
#' @return  a list with the posteriors for the SCR detection function and density covariate parameters (out), s, z
#' @author Ben Augustine
#' @export

mcmc.SCRIPP <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2,beta0=0.05,beta1=0.05),cellArea,Rcpp=TRUE){
    if(Rcpp==TRUE){
      out2=SCRIPPmcmcRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,proppars=proppars,cellArea,keepACs=TRUE)
    }else{
      out2=SCRIPPmcmc(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,proppars=proppars,cellArea,keepACs=TRUE)
    }
    if(keepACs==TRUE){
      list(out=out2$out, sxout=out2$sxout, syout=out2$syout, zout=out2$zout)
    }else{
      list(out=out2$out)
    }
  }
