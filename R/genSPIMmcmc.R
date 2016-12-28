genSPIMmcmc <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=inits,proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE){
    ###
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
    #build possibly true data set
    y.true=y
    for(i in 1:nUnk){
      y.true[ID[i],,]= y.true[ID[i],,]+y.unk[i,,]
    }
    y.true=apply(y.true,c(1,2),sum)
    y.unk=apply(y.unk,c(1,2),sum)

    z=1*(apply(y.true,1,sum)>0)
    z[sample(which(z==0),sum(z==0)/2)]=1 #switch some uncaptured z's to 1.  1/3 is arbitrary. smarter way?
    known.vector=c(rep(1,max(ID)),rep(0,M-max(ID)))

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
      for(l in 1:nUnk){
        #find who you can swap IDs to
        #Should we exlude z=0 guys?
        dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
        possible<- which(dv < swap.tol)+n
        possible=possible[-which(possible==ID[l])]
        possible=possible[z[possible]==1]
        if(length(possible)>0){
          legal=rep(TRUE,length(possible))
          #Check to see if you can swap into these clusters
          for(i in 1:length(possible)){
            check=which(ID==possible[i])#Who else is currently assigned this ID?
            if(length(check)>0){#if false, augmented guy and legal stays true
              for(j in 1:length(check)){
                if(constraints[l,check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                  legal[i]=FALSE
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
            swapped=c(ID[l],newID[l])
            trash<-  sqrt( (s[newID[l],1]- s[(n+1):M,1])^2 + (s[newID[l],2] - s[(n+1):M,2])^2 )
            trash<-   which(dv < swap.tol)+n
            trash=trash[-which(trash==newID[l])]
            if(length(trash)>0){
              legal=rep(TRUE,length(trash))
              #Check to see if you can swap into these clusters
              for(i in 1:length(trash)){
                check=which(ID==trash[i])
                if(length(check)>0){#if false, augmented guy and legal stays true
                  for(j in 1:length(check)){
                    if(constraints[l,check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                      legal[i]=FALSE
                    }
                  }
                }
              }
            }
            trash=trash[legal]
            jump.back<-  1/length(trash)
            #Update y.tmp, z and s for ID[l] and newID
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

            #propose new s
            # Scand=matrix(0,nrow=2,ncol=2)
            # dtmp=matrix(0,nrow=2,ncol=J)
            # Stmp=s[swapped,]
            # lamd.tmp=lamd.tmp.cand=lamd[swapped,]
            # pd.tmp=pd.tmp.cand=pd[swapped,]
            # ll.y.tmp=ll.y.tmp.cand=dbinom(y.cand,K,pd.tmp,log=TRUE)
            # dtmp.cand=dtmp=D[swapped,]
            # for(j in 1:stune){
            #   for(i in 1:2){
            #     inbox=FALSE
            #     while(inbox==FALSE){#propose until cand is in
            #       Scand[i,] <- c(rnorm(1, Stmp[i, 1], proppars$sx2), rnorm(1, Stmp[i, 2], proppars$sy2))
            #       if(useverts==FALSE){
            #         inbox <- Scand[i,1] < xlim[2] & Scand[i,1] > xlim[1] & Scand[i,2] < ylim[2] & Scand[i,2] > ylim[1]
            #       }else{
            #         inbox=inout(Scand[i,],vertices)
            #       }
            #     }
            #     dtmp.cand[i,] <- sqrt((Scand[i,1] - X[, 1])^2 + (Scand[i,2] - X[, 2])^2)
            #     lamd.tmp.cand[i,]<- lam0*exp(-dtmp.cand[i,]*dtmp.cand[i,]/(2*sigma*sigma))
            #     pd.tmp.cand[i,]=1-exp(-lamd.tmp.cand[i,])
            #     ll.y.tmp.cand[i,]=dbinom(y.cand[i,],K,pd.tmp.cand[i,],log=TRUE)
            #     if (runif(1) < exp(sum(ll.y.tmp.cand[i,]) - sum(ll.y.tmp[i,]))) {
            #       Stmp[i,]=Scand[i,]
            #       ll.y.tmp[i,]=ll.y.tmp.cand[i,]
            #       dtmp[i,]=dtmp.cand[i,]
            #       lamd.tmp[i,]=lamd.tmp.cand[i,]
            #       pd.tmp[i,]=pd.tmp.cand[i,]
            #     }
            #   }
            # }
            Scand=matrix(0,nrow=2,ncol=2)
            for(j in 1:stune){
              for(i in 1:2){
                inbox=FALSE
                while(inbox==FALSE){#propose until cand is in
                  Scand[i,] <- c(rnorm(1, Stmp[i, 1], proppars$sx2), rnorm(1, Stmp[i, 2], proppars$sy2))
                  if(useverts==FALSE){
                    inbox <- Scand[i,1] < xlim[2] & Scand[i,1] > xlim[1] & Scand[i,2] < ylim[2] & Scand[i,2] > ylim[1]
                  }else{
                    inbox=inout(Scand[i,],vertices)
                  }
                }
                dtmp.cand[i,] <- sqrt((Scand[i,1] - X[, 1])^2 + (Scand[i,2] - X[, 2])^2)
                lamd.tmp.cand[i,]<- lam0*exp(-dtmp.cand[i,]*dtmp.cand[i,]/(2*sigma*sigma))
                pd.tmp.cand[i,]=1-exp(-lamd.tmp.cand[i,])
                ll.y.tmp.cand[i,]=dbinom(y.cand[i,],K,pd.tmp.cand[i,],log=TRUE)
                if (runif(1) < exp(sum(ll.y.tmp.cand[i,]) - sum(ll.y.tmp[i,]))) {
                  Stmp[i,]=Scand[i,]
                  ll.y.tmp[i,]=ll.y.tmp.cand[i,]
                  dtmp[i,]=dtmp.cand[i,]
                  lamd.tmp[i,]=lamd.tmp.cand[i,]
                  pd.tmp[i,]=pd.tmp.cand[i,]
                }
              }
            }
            if(runif(1)<exp(sum(ll.y.tmp)-sum(ll.y[swapped,]))*(jump.back/jump.probability)){
              y.true[swapped,]=y.cand
              ll.y[swapped,]=ll.y.tmp
              s[swapped, ]=Stmp
              D[swapped, ]=dtmp
              lamd[swapped, ]=lamd.tmp
              pd[swapped,]=pd.tmp
              ID=newID
              known.vector[swapped]=1*(rowSums(y.cand)>0)
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
      list(out=out, sxout=sxout, syout=syout, zout=zout)
    }else{
      list(out=out)
    }
  }

