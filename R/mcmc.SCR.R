#' Run MCMC algorithm for basic SCR model.
#' @param data a list produced by simSCR or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam0=lam0,sigma=sigma)
#' @param proppars a list of tuning parameters for the proposal distributions
#' @param storeLatent a logical indicating whether or not to keep the posteriors for z and s
#' @param Rcpp a logical indicating whether or not to use Rcpp
#' @return  a list with the posteriors for the SCR parameters (out), s, z
#' @author Ben Augustine, Andy Royle
#' @description This function runs the MCMC algorithm for the basic SCR model.  The data list should have the following elements:
#' 1.  y, a n x J capture history
#' 2.  X,  a matrix with the X and Y trap locations in the first two columns
#' 3. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' @export

mcmc.SCR <-
function(data,niter=1000,nburn=0, nthin=1,M = M, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),storeLatent=TRUE,Rcpp=TRUE){
  if(Rcpp==TRUE){ #Do we use Rcpp?
    out2=SCRmcmcRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,proppars=proppars,storeLatent=storeLatent)
  }else{#Don't use Rcpp
    out2=SCRmcmc(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,proppars=proppars,storeLatent=storeLatent)
  }
  if(storeLatent==TRUE){
    list(out=out2$out, sxout=out2$sxout, syout=out2$syout, zout=out2$zout, ID_Lout=out2$ID_Lout,ID_Rout=out2$ID_Rout)
  }else{
    list(out=out2$out)
  }
}

