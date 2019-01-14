#' Spatial partial identity MCMC algorithm.
#' @param data a list produced by sim2side or in the same format
#' @param niter number of MCMC iterations to run
#' @param  nburn number of MCMC iterations to discard as burn in
#' @param nthin MCMC thinning parameter. Record output on every nthin iterations.  nthin=1 corresponds to no thinning
#' @param M The size of the augmented superpopulation
#' @param inits a list of user-supplied initial values.  inits=list(psi=psi,lam01=lam01,lam02=lam02,sigma=sigma)
#' @param  swap number of IDs to swap on each MCMC iteration
#' @param  swap.tol the search radius within which to search for partial ID activity centers to match with
#' @param proppars a list of tuning parameters for the proposal distributions
#' @param storeLatent a logical indicating whether or not to keep the posteriors for z, s, ID_L, and ID_R
#' @param Rcpp a logical indicating whether or not to use Rcpp
#' @return  a list with the posteriors for the SCR parameters (out), s, z, ID_L and ID_R
#' @author Ben Augustine, Andy Royle
#' @description This function runs the MCMC algorithm for the spatial partial identity model.  The data list should have the following elements:
#' 1.  both, a n_both x 3 x K x J both side data array.  If n_both=0 as in an all single camera study, the first dimension is 0 and the data
#' set should still have 4 dimensions, 0 x 3 x K x J.
#' 2.  left, a n_left x 3 x K x J left side data array.
#' 3.  right, a n_right x 3 x K x J left side data array.
#' 4.  IDknown a vector listing the index of complete identity individuals.  It is assumed individuals are sorted such that the complete identity
#' individuals are listed first.  So if there are 7 complete identity individuals, IDknown=1:7.
#' 5. X a matrix with the X and Y trap locations in the first two columns and the number of cameras (1 or 2) at each trap in the third.
#' 6. either buff or vertices.  buff is the fixed buffer for the traps to produce the state space.  It is applied to the minimum and maximum
#' X and Y locations, producing a square or rectangular state space.  vertices is a matrix with the X and Y coordinates of a polygonal state
#' space.
#' 7. an optional tf ,a vector or matrix indicating trap operation. If not accounting for operation across occasions, 
#' tf is a 1 x J vector indicating the number of occasions each trap was operational.  In this scenario,
#' single or double camera stations are either on or off.  If accounting for operation across occasions, tf is a
#' J x K matrix with entries 2 if 2 cameras were operational, 1 if a single camera was operational, and 0 if no
#' cameras were operational.
#' @examples
#' \dontrun{N=50
#'p01=0.13
#'p02=0.2
#'lam01=-log(1-p01)
#'lam02=-log(1-p02)
#'sigma=0.50
#'K=5
#'buff=2
#'niter=1000 #should run more than this and discard a burn in
#'nburn=1
#'nthin=1
#'xlim<- c(1,10)
#'ylim<- c(1,10)
#'X<- expand.grid(3:8,3:8) #6x6 trapping array
#'X=cbind(X,1) #add number of cameras at each trap
#'X[which(X[,2]%in%c(4,7)),3]=2 #switch the second and fifth row of traps to double cameras
#'#Simulate some data
#'data=sim2side(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff)
#'inits=list(psi=0.5,lam01=lam01,lam02=lam02,sigma=sigma)
#'a=Sys.time()
#'store=mcmc.2side(data,niter=niter,nburn=nburn,nthin=nthin, M = 100,inits=inits,swap=10)
#'b=Sys.time()
#'b-a
#'#plot posteriors
#'par(mfrow=c(2,2))
#'plot(store$out[,1],main="lam01",type="l")
#'plot(store$out[,2],main="lam02",type="l")
#'plot(store$out[,3],main="sigma",type="l")
#'plot(store$out[,4],main="N",type="l")
#'par(mfrow=c(1,1))}
#' @export

mcmc.2side <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200, inits=inits,swap=10,swap.tol=1,proppars=list(lam01=0.05,lam02=0.05,sigma=0.1,sx=0.2,sy=0.2),
           storeLatent=TRUE,Rcpp=TRUE){
    if(Rcpp==TRUE){ #Do we use Rcpp?
      if("tf"%in%names(data)){ #Do we have a trap operation file?
        if(length(dim(data$tf))==2){ #Is trap file 2D?
          if(!any(data$tf==2)){
            stop("Collapse 2-D trap file to 1-D vector of counts for days on at each trap because there are no 2's in the current tf")
          }
          out2=mcmc.2sidetfFullRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,swap=swap,swap.tol=swap.tol,
                                    proppars=proppars,storeLatent=storeLatent)
        }else{ #1D trap file
          out2=mcmc.2sideRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,swap=swap,swap.tol=swap.tol,
                            proppars=proppars,storeLatent=storeLatent)
        }
      }else{#No trap file
        out2=mcmc.2sideRcpp(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,swap=swap,swap.tol=swap.tol,
                            proppars=proppars,storeLatent=storeLatent)
      }
    }else{#Don't use Rcpp
      if(length(dim(data$tf)==2)){ #Is trap file 2D?
        out2=mcmc.2sidetf(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,swap=swap,swap.tol=swap.tol,
                          proppars=proppars,storeLatent=storeLatent)
      }else{#trap file 1D
        out2=mcmc.2sideR(data,niter=niter,nburn=nburn, nthin=nthin, M = M, inits=inits,swap=swap,swap.tol=swap.tol,
                          proppars=proppars,storeLatent=storeLatent)
      }
    }
    if(storeLatent==TRUE){
      list(out=out2$out, sxout=out2$sxout, syout=out2$syout, zout=out2$zout, ID_Lout=out2$ID_Lout,ID_Rout=out2$ID_Rout)
    }else{
      list(out=out2$out)
    }
  }
