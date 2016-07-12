e2dist<-function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

cellprobs<- function(lamd){
  # For gaussian hazard model convert lamda(s,x) to p(s,x)
  N<- dim(lamd)[1]
  J<- dim(lamd)[2] # traps
  # just 2 outcomes: left captures and right captures
  pmat<- array(NA,dim=c(N,J,2))
  for(j in 1:J){
    pmat[,j,]<- 1-exp(-lamd[,j,])
  }
  pmat
}

#' Simulate data from camera trap SCR study with partial identity
#' @param N a vector indicating the number of individuals to simulate
#' @param lam01 the single side detection function hazard rate
#' @param lam02 the double side detection function hazard rate
#' @param sigma the spatial scale parameter
#' @param K the number of capture occasions
#' @param X the J x 3 matrix of trap locations and number of cameras at each trap. columns are X, Y, #cams (1 or 2)
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param obstype observation type, either "bernoulli" or "poisson"
#' @param tf an optional J x K matrix of dynamic trap operation indicating when and where 0, 1, or 2 cams were deployed
#' @return a list containing the both, left, and right data sets, the both plus left only and both plus right
#' only data sets (stored in y.obs), the activity centers, the trap object, and several other data objects and summaries.
#' @description This function simulates data from a camera trap SCR study where some or all identities are partial.  The trap
#' matrix, X, should contain the X and Y locations of traps in columns 1 and 2 and the number of cameras at each
#' site in column 3 (either 1 or 2).  The extent of the state space is controlled by "buff", which buffers the
#' minimum and maximum X and Y extents.  Therefore, it probably only makes sense for square or rectangular grids.  Functionality
#' for user-provided polygonal state spaces will be added in the future.
#' @author Ben Augustine
#' @export

sim2side <-
  function(N=50,lam01=0.1,lam02=0.2,sigma=0.50,K=10,X=X,buff=3,obstype="bernoulli",tf=NA){
    library(abind)
    J<- nrow(X)
    #Check for user errors (add more later)
    if(ncol(X)!=3){
      stop("X must have 3 columns, X, Y, # cameras")
    }
    if(!is.null(dim(tf))){
      if(!all(dim(tf)==c(J,K))){
        stop("tf must be of dimension J x K")
      }
      if(!all(tf==1|tf==2|tf==0)){
        stop("elements of tf must be 0, 1, or 2")
      }
      twos=which(apply(tf,1,function(x){any(x==2)}))
      if(!all(twos%in%which(X[,3]==2))){
        stop("trap object X must have a 2 in 3rd column if tf indicates it had 2 cams operational at some point during survey")
      }
      if(any(twos%in%which(X[,3]==1))){
        stop("trap object must have a 1 in 3rd column if tf indicates it never had 2 cams operational")
      }
      usetf=TRUE
    }else{
      usetf=FALSE
      tf=matrix(rep(X[,3],K),ncol=K)
    }
    #######Capture process######################
    # simulate a population of activity centers
    s<- cbind(runif(N, min(X[,1])-buff,max(X[,1])+buff), runif(N,min(X[,2])-buff,max(X[,2])+buff))
    D<- e2dist(s,X)
    lamd<- abind(lam01*exp(-D*D/(2*sigma*sigma)),lam02*exp(-D*D/(2*sigma*sigma)),along=3)

    # Initialize the encounter history data for left and right sides
    left <-right <- both <-array(0,dim=c(N,K,J))
    if(obstype=="bernoulli"){
      pd=cellprobs(lamd)
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            #P(A or B)=P(A)+P(B)-P(A and B)
            left[i,k,j]=rbinom(1,1,(tf[j,k]==1)*pd[i,j,1]+(tf[j,k]==2)*(2*pd[i,j,1]-pd[i,j,1]^2)) #single side p. two chances for capture with 2 cameras
            right[i,k,j]=rbinom(1,1,(tf[j,k]==1)*pd[i,j,1]+(tf[j,k]==2)*(2*pd[i,j,1]-pd[i,j,1]^2)) #single side p
            both[i,k,j]=rbinom(1,1,pd[i,j,2])*(tf[j,k]==2)*1  #both side lambda multiplied by indicator for 2 traps at site
          }
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            left[i,k,j]=rpois(1,tf[j,k]*lamd[i,j,1]) #single side lambda multiplied by number of traps at site
            right[i,k,j]=rpois(1,tf[j,k]*lamd[i,j,1]) #single side lambda multiplied by number of traps at site
            both[i,k,j]=rpois(1,lamd[i,j,2])*(tf[j,k]==2)*1  #both side lambda multiplied by indicator for 2 traps at site
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
    if(usetf==TRUE){
      out<-list(left=left,right=right,both=both,y.obs=y.obs,s=s,X=X,tf=tf,IDknown=IDknown,n=n,nside=nside,Srecap=Srecap,ID_L=ID_L,ID_R=ID_R, K=K,buff=buff,obstype=obstype)
    }else{
      out<-list(left=left,right=right,both=both,y.obs=y.obs,s=s,X=X,IDknown=IDknown,n=n,nside=nside,Srecap=Srecap,ID_L=ID_L,ID_R=ID_R, K=K,buff=buff,obstype=obstype)
    }
    return(out)
  }
