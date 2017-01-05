genSPIMmcmc <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=inits,
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,swap.tol=1,s.lim=6,res=0.2){
    library(abind)
    y<-data$yID
    y.unk<-data$yUnkobs
    X<-as.matrix(data$X)
    J<-nrow(X)
    n<- dim(y)[1]
    K=data$K
    constraints=data$constraints

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
    y<- abind(y,array(0, dim=c( M-dim(y)[1],J, K)), along=1)

    ##match up unknown but constrained samples to true z's
    ##Currently aggregating allowed matches caught at same trap
    if(!any(is.na(constraints))){
      nUnk=nrow(constraints)
      ID=rep(NA,nUnk)
      assign=nUnk
      idx=1:assign
      IDassign=1:assign+n
      used=c()
      while(assign>0){
        focaltrap=which(rowSums(y.unk[idx[1],,])>0)
        cands=which(constraints[idx[1],]==1)
        cands=setdiff(cands,used)
        #cands may be consistent with focal but not each other
        if(length(cands)>1){
          cands=cands[which(rowSums(constraints[cands,cands])==length(cands))]
        }
        keep=rep(0,length(cands))
        for(j in 1:length(cands)){
          candtraps=which(rowSums(y.unk[cands[j],,])>0)
          if(focaltrap%in%candtraps){
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
    }else{
      nUnk=0
    }
    #build possibly true data set
    y.true=y
    if(nUnk>0){
      for(i in 1:nUnk){
        y.true[ID[i],,]= y.true[ID[i],,]+y.unk[i,,]
      }
    }
    y.true=apply(y.true,c(1,2),sum)
    y.unk=apply(y.unk,c(1,2),sum)
    y=apply(y,c(1,2),sum)

    z=1*(apply(y.true,1,sum)>0)
    z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  1/3 is arbitrary. smarter way?
    if(nUnk>0){
      known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))
    }else{
      known.vector=c(rep(1,n),rep(0,M-n))
    }

    #Optimize starting locations given where they are trapped.
    s<- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2])) #assign random locations
    idx=which(rowSums(y.true)>0) #switch for those actually caught
    for(i in idx){
      trps<- X[y.true[i,]>0,1:2]
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
    if(nUnk>0){
      locs=apply(y.unk,1,function(x){which(x>0)})
      locs=X[locs,]
    }

    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=3)
    dimnames(out)<-list(NULL,c("lam0","sigma","N"))
    sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
    IDout=matrix(NA,nrow=nstore,ncol=length(ID))
    idx=1 #for storing output not recorded every iteration

    D<- e2dist(s, X)
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
    pd=pd.cand=1-exp(-lamd)
    ll.y=ll.y.cand= dbinom(y.true,K,pd*z,log=TRUE)

    for(iter in 1:niter){
      #Update lam0
      ll.y= dbinom(y.true,K,pd*z,log=TRUE)
      llysum=sum(ll.y)
      lam0.cand<- rnorm(1,lam0,proppars$lam0)
      if(lam0.cand > 0){
        lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
        pd.cand=1-exp(-lamd.cand)
        ll.y.cand= dbinom(y.true,K,pd.cand*z,log=TRUE)
        llycandsum=sum(ll.y.cand)
        if(runif(1) < exp(llycandsum-llysum)){
          lam0<- lam0.cand
          lamd=lamd.cand
          pd=pd.cand
          ll.y=ll.y.cand
          llysum=llycandsum
        }
      }
      #Update sigma
      sigma.cand<- rnorm(1,sigma,proppars$sigma)
      if(sigma.cand > 0){
        lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
        pd.cand=1-exp(-lamd.cand)
        ll.y.cand= dbinom(y.true,K,pd.cand*z,log=TRUE)
        llycandsum=sum(ll.y.cand)
        if(runif(1) < exp(llycandsum-llysum)){
          sigma<- sigma.cand
          lamd=lamd.cand
          pd=pd.cand
          ll.y=ll.y.cand
          llysum=llycandsum
        }
      }
      #Update ID
      if(nUnk>0){
        for(l in 1:nUnk){
          #find who you can swap IDs to
          #Should we exclude z=0 guys?
          dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
          possible<- which(dv < swap.tol)+n
          possible=possible[-which(possible==ID[l])]
          possible=possible[z[possible]==1]
          if(length(possible)>0){
            thiscluster=which(ID==ID[l])
            legal=rep(TRUE,length(possible))
            #Check to see if you can swap into these clusters
            for(i in 1:length(possible)){
              check=which(ID==possible[i])#Who else is currently assigned this possible new ID?
              if(length(check)>0){#if false, no ID assigned to this guy and legal stays true
                for(k in 1:length(thiscluster)){#loop through guys in focal cluster
                  for(j in 1:length(check)){#check against guys in possible cluster
                    if(constraints[thiscluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                      legal[i]=FALSE
                    }
                  }
                }
              }
            }
            possible=possible[legal]
            if(length(possible)>0){#Can update
              jump.probability<- 1/length(possible) # h(theta*|theta)
              newID=ID
              if(length(possible)>1){
                newID[l] <-  sample( possible, 1)
              }else{
                newID[l]=possible
              }
              trash<-  sqrt( (s[newID[l],1]- s[(n+1):M,1])^2 + (s[newID[l],2] - s[(n+1):M,2])^2 )
              trash<-   which(dv < swap.tol)+n
              trash=trash[-which(trash==newID[l])]
              legal=rep(TRUE,length(trash))
              backcluster=which(newID==newID[l])#which guys are in the new ID cluster?
              #Check to see if you can swap into these clusters
              for(i in 1:length(trash)){
                check=which(newID==trash[i])
                for(k in 1:length(backcluster)){#loop through guys in back cluster
                  if(length(check)>0){#if false, no ID assigned to this back guy and legal stays true
                    for(j in 1:length(check)){#check against back guys in proposed cluster
                      if(constraints[backcluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                        legal[i]=FALSE
                      }
                    }
                  }
                }
              }
              trash=trash[legal]
              jump.back<-  1/length(trash)
              swapped=c(ID[l],newID[l])#order swap.out then swap.in
              cluster1guys=which(newID==swapped[1])
              cluster2guys=which(newID==swapped[2])
              y.cand=matrix(0,nrow=2,ncol=J)
              if(length(cluster1guys)>0){
                if(length(cluster1guys)>1){
                  y.cand[1,]=colSums(y.unk[cluster1guys,])
                }else{
                  y.cand[1,]=y.unk[cluster1guys,]
                }
              }
              if(length(cluster2guys)>0){
                if(length(cluster2guys)>1){
                  y.cand[2,]=colSums(y.unk[cluster2guys,])
                }else{
                  y.cand[2,]=y.unk[cluster2guys,]
                }
              }
              #Propose new S from full conditional
              s.cand=matrix(NA,nrow=2,ncol=2)
              dtmp=matrix(NA,nrow=2,ncol=J)
              for(i in 1:2){
                grid.locs=as.matrix(expand.grid(seq(s[swapped[i],1]-s.lim/2,s[swapped[i],1]+s.lim/2,res),
                                 seq(s[swapped[i],2]-s.lim/2,s[swapped[i],2]+s.lim/2,res)))
                #Check for any off state space
                rem=which(grid.locs[,1]<xlim[1]|grid.locs[,1]>xlim[2]|grid.locs[,2]<ylim[1]|grid.locs[,2]>ylim[2])
                if(length(rem)>0){
                  grid.locs=grid.locs[-rem,]
                }
                dtmp.tmp=e2dist(grid.locs,X)
                lamd.tmp=lam0*exp(-dtmp.tmp*dtmp.tmp/(2*sigma*sigma))
                pd.tmp=1-exp(-lamd.tmp)
                npoints=nrow(grid.locs)
                grid.lik=rep(NA,nrow=npoints)
                for(j in 1:npoints){
                  grid.lik[j]=sum(dbinom(y.cand[i,],K,pd.tmp[j,]))
                }
                grid.probs=grid.lik/sum(grid.lik)
                choose=sample(1:npoints,1,prob=grid.probs)
                s.cand[i,]=grid.locs[choose,]
                dtmp[i,]=dtmp.tmp[i,]
                lamd.cand[swapped[i],]=lamd.tmp[choose,]
                pd.cand[swapped[i],]=pd.tmp[choose,]
              }
# image(seq(s[swapped[i],1]-s.lim/2,s[swapped[i],1]+s.lim/2,res),
#       seq(s[swapped[i],2]-s.lim/2,s[swapped[i],2]+s.lim/2,res),matrix(grid.probs,nrow=sqrt(npoints)))


              ll.y.cand[swapped,]=dbinom(y.cand,K,pd.cand[swapped,],log=TRUE)
              #Check for errors
              for(j in 1:length(newID)){
                focal=ID[j]
                thiscluster=which(newID==focal)
                if(length(thiscluster)>1){
                  for(k in 1:length(thiscluster))
                    if(any(constraints[thiscluster,thiscluster]==0)){
                      stop("error")
                    }
                }
              }
              if(runif(1)<exp(sum(ll.y.cand[swapped,])-sum(ll.y[swapped,]))*(jump.back/jump.probability)){
                y.true[swapped,]=y.cand
                D[swapped,]=dtmp
                lamd[swapped,]=lamd.cand[swapped,]
                pd[swapped,]=pd.cand[swapped,]
                s[swapped,]=s.cand
                ll.y[swapped,]=ll.y.cand[swapped,]
                ID=newID
                known.vector[swapped]=1*(rowSums(y.cand)>0)
              }
            }
          }
        }
      }
      #Update psi gibbs
      ## probability of not being captured in a trap AT ALL
      pbar=(1-pd)^K
      prob0<- exp(rowSums(log(pbar)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      ll.y= dbinom(y.true,K,pd*z,log=TRUE)
      psi <- rbeta(1, 1 + sum(z), 1 + M - sum(z))

      ## Now we have to update the activity centers
      for (i in 1:M) {
        Scand <- c(rnorm(1, s[i, 1], proppars$sx), rnorm(1, s[i, 2], proppars$sy))
        if(useverts==FALSE){
          inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
        }else{
          inbox=inout(Scand,vertices)
        }
        if (inbox) {
          dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
          lamd.cand[i,]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
          pd.cand[i,]=1-exp(-lamd.cand[i,])
          ll.y.cand[i,]= dbinom(y.true[i,],K,pd.cand[i,]*z[i],log=TRUE)
          if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
            s[i, ]=Scand
            D[i, ]=dtmp
            lamd[i, ]=lamd.cand[i,]
            pd[i,]=pd.cand[i,]
            ll.y[i,]=ll.y.cand[i,]
          }
        }
      }
      #Do we record output on this iteration?
      if(i>nburn&i%%nthin==0){
        sxout[idx,]<- s[,1]
        syout[idx,]<- s[,2]
        zout[idx,]<- z
        IDout[idx,]=ID
        out[idx,]<- c(lam0,sigma ,sum(z))
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else{
      list(out=out)
    }
  }

