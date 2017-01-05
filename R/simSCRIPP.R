#' Simulate data from camera trap SCR study with inhomogenous density in discrete space and one categorical covariate with 2 levels
#' @param N a vector indicating the number of individuals to simulate
#' @param lam0 the detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the K x 2 matrix of trap locations
#' @param grid the data frame holding the state space points and covariate values. Should have columns labeled x, y, and cov with cov==1 or 2
#' @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
#' @description This function simulates data from a camera trap SCR study with inhomogenous density for two categorical covariates.
#' Binomial observation model.  Will generalize later.
#' @author Ben Augustine
#' @export

simSCRIPP <- function(N=120,lam0=0.2,sigma=0.50,K=10,X=X,grid,Dparms,plot=TRUE){
    ############Set up density covariate
    beta=Dparms[1]
    grid$cp=exp((0+beta*(grid$cov==2))) / sum(0+exp(beta*(grid$cov==2)))
    #s.tmp=rmultinom(1, N, grid$cp) # a single realization to be ignored
    npix=nrow(grid)

    # # simulate a population of activity centers
    s=matrix(NA,nrow=N,ncol=3)
    for(i in 1:N) {
      s.i <- sample(1:npix, 1, prob=grid$cp)
      sx <- grid[s.i, "x"]
      sy <- grid[s.i, "y"]
      s[i,] <- c(s.i, sx, sy)
    }
    if(plot==TRUE){
      grid2=xtabs(grid$cov~grid$x+grid$y)
      image(x=as.numeric(rownames(grid2)),y=as.numeric(colnames(grid2)),z=grid2,xlab="X",ylab="Y")
      points(s[,2:3])
      points(X,pch=4)
    }
    D<- e2dist(s[,2:3],X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)

    #######Capture process######################
    # Simulate encounter history
    y <-array(0,dim=c(N,K,J))
    pd=cellprobsSCR(lamd)
    for(i in 1:N){
      for(j in 1:J){
        for(k in 1:K){
          y[i,k,j]=rbinom(1,1,pd[i,j])
        }
      }
    }
    y=y[-which(rowSums(y)==0),,]
    n=dim(y)[1]
    #Count spatial recaps
    y2D=apply(y,c(1,3),sum)
    scaps=rowSums(1*(y2D>0))
    scaps[scaps>0]=scaps[scaps>0]-1
    nscap=sum(scaps>0)
    sumscap=sum(scaps)
    out<-list(y=y,s=s,X=X, K=K,n=n,nscap=nscap,sumscap=sumscap,Dparms=Dparms,grid=grid)
    return(out)
  }
