#' Run MCMC algorithm for the naive independence estimator.
#' @param data a list produced by sim2side or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam01=lam01,sigma=sigma)
#' @param proppars a list of tuning parameters for the proposal distributions
#' @param keepACs a logical indicating whether or not to return the posteriors for z, and s
#' @return  a list with the posteriors for the SCR parameters (out), s, z
#' @author Ben Augustine, Andy Royle
#' @description This function fits the naive independence estimator from Augustine et al. 2016.  The data list should have the following elements:
#' 1.  both, a n_both x 3 x K x J left side data array.  Not required if there are no double camera stations.
#' 2.  left, a n_left x 3 x K x J left side data array.
#' 3.  right, a n_right x 3 x K x J left side data array.
#' 4. X a matrix with the X and Y trap locations in the first two columns and the number of cameras (1 or 2) at each trap in the third.
#' 5. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' inits needs elements for psi1, psi2, and psi3 if there are some double camera stations or psi1 and psi2 if there are only single camera stations.
#' proppars needs elements for lam01 and lam02 if there are some double camera stations or lam01 if there are only single camera stations.
#' @export

mcmc.2side.ind <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,Rcpp=TRUE){
    if(Rcpp==TRUE){#Use Rcpp?
      if(any(data$X==2)){#Are there any double cam stations?
        if(!all(c("psi1","psi2","psi3")%in%names(inits))){
          stop("Initial values for psi1, psi2, and psi3 need to be specified")
        }
        out2=mcmc.2side.ind3Rcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,
                                 proppars=proppars,keepACs=keepACs)
      }else{#all single cams
        if(!all(c("psi1","psi2")%in%names(inits))){
          stop("Initial values for psi1 and psi2 need to be specified")
        }
        if("lam02"%in%names(proppars)){
          proppars[-which(names(proppars)=="lam02")]
        }
        out2=mcmc.2side.ind2Rcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,
                                 proppars=proppars,keepACs=keepACs)
      }
    }else{#No Rcpp
      if(any(data$X==2)){#Are there any double cam stations?
        if(!all(c("psi1","psi2","psi3")%in%names(inits))){
          stop("Initial values for psi1, psi2, and psi3 need to be specified")
        }
        out2=mcmc.2side.ind3(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,
                                 proppars=proppars,keepACs=keepACs)
      }else{#all single cams
        if(!all(c("psi1","psi2")%in%names(inits))){
          stop("Initial values for psi1 and psi2 need to be specified")
        }
        if("lam02"%in%names(proppars)){
          proppars[-which(names(proppars)=="lam02")]
        }
        out2=mcmc.2side.ind2(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,
                                 proppars=proppars,keepACs=keepACs)
      }
    }
    if(keepACs==TRUE){
      list(out=out2$out, sLxout=out2$sLxout, sLyout=out2$sLyout,sRxout=out2$sRxout, sRyout=out2$sRyout, zLout=out2$zLout,zRout=out2$zRout)
    }else{
      list(out=out2$out)
    }
  }
