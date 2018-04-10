genSPIMmcmcRcpp <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=inits,obstype="bernoulli",
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,swap.tol=1){
    library(abind)
    y<-data$yID
    y.unk<-data$yUnkobs
    X<-as.matrix(data$X)
    J<-nrow(X)
    n<- dim(y)[1]
    K=data$K
    constraints=data$constraints
    
    #data checks
    if(length(dim(y))!=3){
      stop("dim(yID) must be 3 even if no guys captured")
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
    ##pull out initial values
    psi<- inits$psi
    lam0<- inits$lam0
    sigma<- inits$sigma
    
    #Augment data and make initial complete data set
    if(length(dim(y))!=3){
      y<- abind(y,array(0, dim=c( M-dim(y)[1],J, K)), along=1)
    }else if(dim(y)[1]>0){
      y<- abind(y,array(0, dim=c( M-dim(y)[1],J, K)), along=1)
    }else{
      y<- array(0,dim=c(M,J,K))
    }
    
    ##match up unknown but constrained samples to true z's
    ##Currently aggregating allowed matches caught at at least 1 same trap
    if(!any(is.na(constraints))){
      nUnk=nrow(constraints)
      ID=rep(NA,nUnk)
      assign=nUnk
      idx=1:assign
      IDassign=1:assign+n
      used=c()
      while(assign>0){
        focaltraps=which(rowSums(y.unk[idx[1],,])>0)
        cands=which(constraints[idx[1],]==1)
        cands=setdiff(cands,used)
        #cands may be consistent with focal but not each other
        if(length(cands)>1){
          cands=cands[which(rowSums(constraints[cands,cands])==length(cands))]
        }
        keep=rep(0,length(cands))
        for(j in 1:length(cands)){
          candtraps=which(rowSums(y.unk[cands[j],,])>0)
          if(any(focaltraps%in%candtraps)){
            keep[j]=1
          }
        }
        cluster=cands[keep==1]
        ID[cluster]=IDassign[1]
        IDassign=IDassign[-1]
        idx=idx[-which(idx%in%cluster)]
        used=c(used,cluster)
        assign=assign-length(cluster)
      }
      if(obstype=="bernoulli"){#break apart clusters that create y_ijk>1
        y.true=y
        if(nUnk>0){
          for(i in 1:nUnk){
            if(ID[i]>M){
              stop("Need to raise M to initialize latent observations")
            }
            y.true[ID[i],,]= y.true[ID[i],,]+y.unk[i,,]
          }
        }
        fix=which(y.true>1,arr.ind=TRUE)
        nfix=nrow(fix)
        while(nfix>0){
          if(is.matrix(fix)==FALSE){
            fix=matrix(fix,ncol=3)
          }
          if(nrow(fix)>1){#multiple rows
            who=which(y.unk[,fix[1,2],fix[1,3]]==1)
            nreassign=length(who)-1
            ID[who[2:(2+nreassign-1)]]=IDassign[1:(1+nreassign-1)]
            IDassign=IDassign[-(1:(1+nreassign-1))]
            fix=fix[-1,]
          }else{
            who=which(y.unk[,fix[2],fix[3]]==1)
            
            nreassign=length(who)-1
            ID[who[2:(2+nreassign-1)]]=IDassign[1:(1+nreassign-1)]
            IDassign=IDassign[-(1:(1+nreassign-1))]
            nfix=0
          }
        }
      }
    }else{
      nUnk=0
      ID=NA
      constraints=matrix(c(1,1,1,1),nrow=2)#gotta feed rcpp a matrix
    }
    #build possibly true data set
    y.true=y
    if(nUnk>0){
      for(i in 1:nUnk){
        if(ID[i]>M){
          stop("Need to raise M to initialize latent observations")
        }
        y.true[ID[i],,]= y.true[ID[i],,]+y.unk[i,,]
      }
    }
    y.true2D=apply(y.true,c(1,2),sum)
    y.unk2D=apply(y.unk,c(1,2),sum)
    y2D=apply(y,c(1,2),sum)

    z=1*(apply(y.true,1,sum)>0)
    z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  1/3 is arbitrary. smarter way?
    if(nUnk>0){
      known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))
    }else{
      known.vector=c(rep(1,n),rep(0,M-n))
    }

    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.true2D)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[y.true2D[i,]>0,1:2]
      trps<-matrix(trps,ncol=2,byrow=FALSE)
      s[i,]<- c(mean(trps[,1]),mean(trps[,2]))
    }
    #check to make sure everyone is in polygon
    if("vertices"%in%names(data)){
      vertices=data$vertices
      useverts=TRUE
    }else{
      useverts=FALSE
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
    D<- e2dist(s, X)
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))

    #calcluate ll.unk
    # njk=apply(y.unk,c(2,3),sum)
    # ll.unk=ll.unk.cand=matrix(0,J,K)
    # unk.guys=(n+1):M
    ID2=ID-1#bc Rcpp starts counting at 0
    
    if(obstype=="bernoulli"){
      store=mcmc_genSPIMbern(lam0, sigma,D,lamd,y.true,y.unk,X,z, s, useverts, vertices, xlim,ylim, known.vector, 
                             psi, proppars$lam0,proppars$sigma, proppars$sx, proppars$sy, niter,nburn,nthin,
                             n,nUnk,swap.tol,constraints,ID2)
    }else{
      store=mcmc_genSPIMpois(lam0, sigma,D,lamd,y.true,y.unk,X,z, s, useverts, vertices, xlim,ylim, known.vector, 
                             psi, proppars$lam0,proppars$sigma, proppars$sx, proppars$sy, niter,nburn,nthin,
                             n,nUnk,swap.tol,constraints,ID2)
    }
    out=store[[1]]
    sxout=store[[2]]
    syout=store[[3]]
    zout=store[[4]]
    IDout=store[[5]]
    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else{
      list(out=out)
    }
  }
