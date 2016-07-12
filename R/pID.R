#' Produce the posteriors for P(ID_L_i = ID_R_j) from the posteriors of ID_L and ID_R
#' @param data a list produced by sim2side or in the same format
#' @param ID_Lout the M x niter posterior of ID_L
#' @param ID_Rout the M x niter posterior of ID_R
#' @param swapped TRUE if mcmc2side swapped L and R sides before MCMC sampling, FALSE otherwise
#' @author Ben Augustine
#' @description This function produces the posterior that left individual X is right individual Y.  For complete identity individuals, left individual X =
#' right individual X with probability 1.
#' @return a list with elements pIDL and pIDR, which are also lists.  Element X of pIDL corresponds to left individual X and contains the right
#' individuals it was matched with on at least one MCMC iteration and the posterior for each match.  Similarly for pIDR.
#' @examples
#' \dontrun{
#' N=50
#'p01=0.13
#'p02=0.2
#'lam01=-log(1-p01)
#'lam02=-log(1-p02)
#'sigma=0.50
#'K=5
#'buff=2
#'xlim<- c(1,10)
#'ylim<- c(1,10)
#'X<- expand.grid(3:8,3:8) #6x6 trapping array
#'X=cbind(X,1) #add number of cameras at each trap
#'X[which(X[,2]%in%c(4,7)),3]=2 #switch the second and fifth row of traps to double cameras
#'#Simulate some data
#'data=sim2side(N=N,lam01=lam01,lam02=lam02,sigma=sigma,K=K,X=X,buff=buff)
#'inits=list(psi=0.5,lam01=lam01,lam02=lam02,sigma=sigma)
#'store=mcmc.2sideRcpp(data,niter=niter,nburn=nburn,nthin=nthin, M = 100,inits=inits,swap=10,keepACs=TRUE)
#'b=Sys.time()
#'pID(data,store$ID_L,store$ID_R)
#'}
pID=function(data,ID_Lout,ID_Rout,swapped=FALSE){
  #Did the mcmc2side function swap sides?
  if(swapped==TRUE){
    store=ID_Lout
    ID_Lout=ID_Rout
    ID_Rout=store
  }
  ID_L=data$ID_L
  ID_R=data$ID_R
  pIDL=vector("list",length(ID_L))
  pIDR=vector("list",length(ID_R))
  niter=nrow(ID_Lout)
  #Do left guys first
  for(i in 1:length(ID_L)){
    match=rep(NA,niter)
    for(j in 1:niter){
      bguy=ID_Lout[j,i] #which both guy is this left guy?
      match[j]=which(ID_Rout[j,]==bguy) #which right guy is assigned to that both guy?
    }
    tab=table(match)
    rID=as.numeric(names(tab))
    pIDL[[i]]=cbind(rID,as.numeric(tab/niter))
    colnames(pIDL[[i]])=c("rID","Prob")
  }
  #Then right guys
  for(i in 1:length(ID_R)){
    match=rep(NA,niter)
    for(j in 1:niter){
      bguy=ID_Rout[j,i] #which both guy is this right guy?
      lguy=which(ID_Lout[j,]==bguy) #which left guy is assigned to that both guy?
      match[j]=lguy
    }
    tab=table(match)
    lID=as.numeric(names(tab))
    pIDR[[i]]=cbind(lID,as.numeric(tab/niter))
    colnames(pIDR[[i]])=c("lID","Prob")
  }
  return(list(pIDL=pIDL,pIDR=pIDR))
}
