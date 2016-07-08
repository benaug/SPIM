# Simulate data from camera trap SCR study with inhomogenous density in discrete space and one categorical covariate with 2 levels
# @param N a vector indicating the number of individuals to simulate
# @param lam0 the detection function hazard rate
# @param sigma the spatial scale parameter
# @param K the number of capture occasions
# @param X the K x 2 matrix of trap locations
# @param grid the data frame holding the state space points and covariate values. Should have columns labeled x, y, and cov with cov==1 or 2
# @return a list containing the capture history, activity centers, trap object, and several other data objects and summaries.
# @description This function simulates data from a camera trap SCR study with inhomogenous density for two categorical covariates.
# Binomial observation model.  Will generalize later.
# @author Ben Augustine
# @export

sim2sideSCRIPP <- function(N=120,lam01=0.1,lam02=0.2,sigma=0.50,K=10,X=X,grid,Dparms,plot=TRUE,obstype="bernoulli"){
    library(abind)
    #Check for user errors (add more later)
    if(ncol(X)!=3){
      stop("X must have 3 columns, X, Y, # cameras")
    }
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
    lamd<- abind(lam01*exp(-D*D/(2*sigma*sigma)),lam02*exp(-D*D/(2*sigma*sigma)),along=3)
    J<- nrow(X)

    #######Capture process######################
    # Initialize the encounter history data for left and right sides
    left <-right <- both <-array(0,dim=c(N,K,J))
    if(obstype=="bernoulli"){
      pd=cellprobs(lamd)
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            left[i,k,j]=rbinom(1,1,(X[j,3]==1)*pd[i,j,1]+(X[j,3]==2)*(2*pd[i,j,1]-pd[i,j,1]^2)) #single side p. two chances for capture with 2 cameras
            right[i,k,j]=rbinom(1,1,(X[j,3]==1)*pd[i,j,1]+(X[j,3]==2)*(2*pd[i,j,1]-pd[i,j,1]^2)) #single side p
            both[i,k,j]=rbinom(1,1,pd[i,j,2])*(X[j,3]==2)*1  #both side lambda multiplied by indicator for 2 traps at site
          }
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            left[i,k,j]=rpois(1,X[j,3]*lamd[i,j,1]) #single side lambda multiplied by number of traps at site
            right[i,k,j]=rpois(1,X[j,3]*lamd[i,j,1]) #single side lambda multiplied by number of traps at site
            both[i,k,j]=rpois(1,lamd[i,j,2])*(X[j,3]==2)*1  #both side lambda multiplied by indicator for 2 traps at site
          }
        }
      }
    }else{
      stop("observation model not recognized")
    }

    ######Process data#############
    IDknown=which(apply(both,1,sum)>0)
    Nknown=length(IDknown)
    n=sum(apply(both+left+right,1,sum)>0)
    nside=c(sum(apply(both,1,sum)>0),sum(apply(left,1,sum)>0),sum(apply(right,1,sum)>0))
    #Count true spatial recaps
    y2D=apply(both+left+right,c(1,3),sum)
    scaps=rowSums(1*(y2D>0))
    scaps[scaps>0]=scaps[scaps>0]-1
    nscap=sum(scaps>0)
    sumscap=sum(scaps)
    Srecap=data.frame(nscap=nscap,sumscap=sumscap)
    # relabel data sets if needed so that the left side always has more individuals than right after removing known IDs
    nl.tmp<- sum(apply(left,1,sum)[apply(both,c(1),sum)==0]>0)
    nr.tmp<- sum(apply(right,1,sum)[apply(both,c(1),sum)==0]>0)
    if(nr.tmp > nl.tmp){
      a<- left
      left<- right
      right<- a
    }
    # Move known individuals(after photo ID) to beginning
    if(Nknown>0&Nknown<N){
      IDunknown=which(is.na(match(1:N,IDknown)))
      newleft<- abind(left[IDknown,,],left[IDunknown,,],along=1)
      newright<- abind(right[IDknown,,],right[IDunknown,,],along=1)
      newboth<- abind(both[IDknown,,],both[IDunknown,,],along=1)
      s=s[c(IDknown,IDunknown),]
      left<- newleft
      right<- newright
      both<- newboth
      #Recalculate after sorting
      IDknown=1:length(IDknown)
      # Sort by number of captures in both data set.
      #These are now defined to be in order by true IDs 1:N
      tb<- apply(both,c(1,3),sum)
      cap.both<- apply(tb,1,sum)[(Nknown+1):N]
      oo1<- c(1:Nknown, Nknown + rev(order(cap.both)) )
      both<- both[oo1,,]
      left<- left[oo1,,]
      right<- right[oo1,,]
      s<- s[oo1,]
      #Sort by number of captures in observed (R) left and right data sets.
      #Not sure if this is necessary, inherited from Andy
      tl<- apply(left,c(1,3),sum)
      tr<- apply(right,c(1,3),sum)
      cap.left<- apply(tl,1,sum)[(Nknown+1):N]
      cap.right<- apply(tr,1,sum)[(Nknown+1):N]
      oo2<- c(1:Nknown, Nknown+ rev(order(cap.left)) )
      oo3<- c(1:Nknown, Nknown + rev(order(cap.right)) )
      left<- left[oo2,,]
      right<- right[oo3,,]
      #Record order of true IDs
      ID_L=(1:N)[oo2]
      ID_R=(1:N)[oo3]
    }else{
      #Record order of true IDs
      ID_L=(1:N)
      ID_R=(1:N)
    }
    #Add all zero dimensions so both left and right can be added together
    add=array(0,dim=c(N,K,J))
    both=abind(both,add,add,along=0)
    left=abind(add,left,add,along=0)
    right=abind(add,add,right,along=0)
    both=aperm(both,c(2,1,3,4))
    left=aperm(left,c(2,1,3,4))
    right=aperm(right,c(2,1,3,4))

    #Remove all 0 capture histories.  Don't remove right and left all 0 histories for known IDs
    both<- both[IDknown,,,]
    if(length(dim(both))==3){#if there is 1 both guy need to keep 4d
      both=array(both,dim=c(1,dim(both)))
    }
    lcap=which(apply(left,1,sum)>0)
    rcap=which(apply(right,1,sum)>0)
    left<- left[sort(unique(c(lcap,IDknown))),,,]
    right<- right[sort(unique(c(rcap,IDknown))),,,]
    if(length(dim(left))==3){#if there is 1 both guy need to keep 4d
      left=array(left,dim=c(1,dim(left)))
    }
    if(length(dim(right))==3){#if there is 1 both guy need to keep 4d
      right=array(right,dim=c(1,dim(right)))
    }
    lorder=sort(unique(c(lcap,IDknown)))
    rorder=sort(unique(c(rcap,IDknown)))

    #Update true ID order for left and right for captured guys
    ID_L=ID_L[lorder]
    ID_R=ID_R[rorder]

    #Build all possible observed data sets and count observed spatial recaps
    B2D=apply(both,c(1,4),sum)
    L2D=apply(left,c(1,4),sum)
    R2D=apply(right,c(1,4),sum)
    known=dim(both)[1]
    if(known>0){
      BLR=1*((both[,1,,]+left[1:known,2,,]+right[1:known,3,,])>0) #keep boths and lefts and rights for boths
      if(length(dim(BLR))==2){
        BLR=array(BLR,dim=c(1,dim(BLR)))
      }
      if(dim(left)[1]>known){
        BLRL=abind(BLR,left[(known+1):dim(left)[1],2,,],along=1) #Add unknown lefts
      }else{
        BLRL=BLR
      }
      if(dim(right)[1]>known){
        BLRR=abind(BLR,right[(known+1):dim(right)[1],3,,],along=1) #Add unknown rights
      }else{
        BLRR=BLR
      }
    }else{ #no known guys
      BLR=both[,1,,]
      if(dim(left)[1]>known){
        BLRL=left[,2,,] #Add unknown lefts
      }else{
        BLRL=BLR
      }
      if(dim(right)[1]>known){
        BLRR=right[,3,,] #Add unknown rights
      }else{
        BLRR=BLR
      }
      if(length(dim(BLRL))==2){
        BLRL=array(BLRL,dim=c(1,dim(BLRL)))
      }
      if(length(dim(BLRR))==2){
        BLRR=array(BLRR,dim=c(1,dim(BLRR)))
      }
    }
    #BLR2D=apply(BLR,c(1,3),sum)
    BLRL2D=apply(BLRL,c(1,3),sum)
    BLRR2D=apply(BLRR,c(1,3),sum)
    if(dim(B2D)[1]>0){
      Bscaps=rowSums(1*(B2D>0))
      Bscaps[Bscaps>0]=Bscaps[Bscaps>0]-1
      Bnscap=sum(Bscaps>0)
      Bsumscap=sum(Bscaps)
    }else{
      Bnscap=0
      Bsumscap=0
    }
    if(dim(L2D)[1]>0){
      Lscaps=rowSums(1*(L2D>0))
      Lscaps[Lscaps>0]=Lscaps[Lscaps>0]-1
      Lnscap=sum(Lscaps>0)
      Lsumscap=sum(Lscaps)
    }else{
      Lnscap=0
      Lsumscap=0
    }
    if(dim(R2D)[1]>0){
      Rscaps=rowSums(1*(R2D>0))
      Rscaps[Rscaps>0]=Rscaps[Rscaps>0]-1
      Rnscap=sum(Rscaps>0)
      Rsumscap=sum(Rscaps)
    }else{
      Rnscap=0
      Rsumscap=0
    }
    if(dim(BLRL2D)[1]>0){
      BLRLscaps=rowSums(1*(BLRL2D>0))
      BLRLscaps[BLRLscaps>0]=BLRLscaps[BLRLscaps>0]-1
      BLRLnscap=sum(BLRLscaps>0)
      BLRLsumscap=sum(BLRLscaps)
    }else{
      BLRLnscap=0
      BLRLsumscap=0
    }
    if(dim(BLRR2D)[1]>0){
      BLRRscaps=rowSums(1*(BLRR2D>0))
      BLRRscaps[BLRRscaps>0]=BLRRscaps[BLRRscaps>0]-1
      BLRRnscap=sum(BLRRscaps>0)
      BLRRsumscap=sum(BLRRscaps)
    }else{
      BLRRnscap=0
      BLRRsumscap=0
    }
    Srecap$Bnscap=Bnscap
    Srecap$Bsumscap=Bsumscap
    Srecap$Lnscap=Lnscap
    Srecap$Lsumscap=Lsumscap
    Srecap$Rnscap=Rnscap
    Srecap$Rsumscap=Rsumscap
    Srecap$BLRLnscap=BLRLnscap
    Srecap$BLRLsumscap=BLRLsumscap
    Srecap$BLRRnscap=BLRRnscap
    Srecap$BLRRsumscap=BLRRsumscap
    dim(BLRR)
    dim(BLRL)
    y.obs=list(BLRL=BLRL,BLRR=BLRR)
    out<-list(left=left,right=right,both=both,y.obs=y.obs,s=s,X=X,IDknown=IDknown,n=n,nside=nside,Srecap=Srecap,ID_L=ID_L,ID_R=ID_R,K=K,grid=grid,Dparms=Dparms,obstype=obstype)
    return(out)
  }
