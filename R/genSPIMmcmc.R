genSPIMmcmc <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=NA,obstype="bernoulli",
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,swap.tol=1,nswap=10,IDupdate="Gibbs"){
    library(abind)
    y<-data$yID
    y.unk<-data$yUnkobs
    X<-as.matrix(data$X)  
    J<-nrow(X)
    n<- dim(y)[1]
    K=data$K
    constraints=data$constraints
    if(obstype=="bernoulli"&(dim(y.unk)[1]>0)){#add constraints that create y_ijk>1
      caps=which(y.unk==1,arr.ind=TRUE)
      for(i in 1:nrow(caps)){
        for(j in 1:nrow(caps)){
          if(i==j)break
          if(all(caps[i,2:3]==caps[j,2:3])){
            constraints[i,j]=0
            constraints[j,i]=0
          }
        }
      }
    }
    #data checks
    if(length(dim(y))!=3){
      stop("dim(yID) must be 3 even if no guys captured")
    }
    if(IDupdate=="MH"){
      if(nswap>nrow(y.unk)){
        warning("nswap is larger than number of partial ID guys. Setting nswap=nrow(y.unk)")
        nswap=nrow(y.unk)
      }
    }
    if(dim(data$yUnk)[1]>0){
      if(nrow(y.unk)!=nrow(constraints)){
        stop("data$yUnkobs must have same number of rows as data$constraints")
      }
      if(nrow(constraints)!=ncol(constraints)){
        stop("identity constraint matrix needs to be symmetric")
      }
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
    theta=inits$theta

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
      ID=1:n
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
    
    # some objects to hold the MCMC simulation output
    nstore=(niter-nburn)/nthin
    if(nburn%%nthin!=0){
      nstore=nstore+1
    }
    out<-matrix(NA,nrow=nstore,ncol=5)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","clusters","theta"))
    sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
    IDout=matrix(NA,nrow=nstore,ncol=length(ID))
    idx=1 #for storing output not recorded every iteration
    
    #collapse data to 2D
    y.true=apply(y.true,c(1,2),sum)
    y.unk=apply(y.unk,c(1,2),sum)
    D<- e2dist(s, X)
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
    
    #more ID initialization
    rowcounts=rowSums(y.true[(n+1):M,])
    allcounts=sum(rowcounts)
    rowsumprobs=(rowSums(lamd[(n+1):M,])/sum(lamd[(n+1):M,]))
    expectedcounts=round(rowsumprobs*allcounts)
    for(swap in 1:15){
      fix=which(rowcounts>expectedcounts)+n
      recruit=which(rowcounts<expectedcounts)+n
      if(length(fix)==0|length(recruit)==0) break  #stop if we run out of guys to fix or recruit
      for(i in fix){
        cands=which(ID==i)#find samples assigned to this guy
        if(length(cands)>1){
          pick=sample(cands,1)#pick one to give away
        }else{
          pick=cands
        }
        dv= sqrt((s[i,1] - s[recruit, 1])^2 + (s[i,2] - s[recruit, 2])^2)
        pick2=which(dv==min(dv))[1] #pick the closest guy to swap ID with
        ID[pick]=recruit[pick2]
        #update data after ID swap
        y.true[i,]=y.true[i,]-y.unk[pick,]
        y.true[recruit[pick2],]=y.true[recruit[pick2],]+y.unk[pick,]
        rowcounts[i-n]=sum(y.true[i,])
        rowcounts[recruit[pick2]-n]=sum(y.true[recruit[pick2],])
        z[recruit[pick2]]=1
      }
    }
    if(max(expectedcounts-rowcounts)>5){
      warning("one guy initialized with >5 counts than expected. Make sure it's being updated.")
    }
    known.vector=1*(rowSums(y.true)>0)
    
    
    if(obstype=="bernoulli"){
      pd=pd.cand=1-exp(-lamd)
      ll.y=ll.y.cand= dbinom(y.true,K,pd*z,log=TRUE)
    }else{
      ll.y=ll.y.cand= dpois(y.true,K*lamd*z,log=TRUE)
    }
    if(nUnk>0&IDupdate=="Gibbs"){
      ntraps=rowSums(y.unk>0)
      if(any(ntraps>1)){#at least one connection across traps
        stop("Can't use Gibbs with connections across traps :(")
       
      }else if(!all(constraints==1)){#no connections across traps
        Gibbs="simple"
      }else{
        Gibbs="RC"
        njGibbs=colSums(y.unk)
        trapcaps=which(njGibbs>0)
      }
    }else if(IDupdate=="MH"){
    }else{
      stop("IDupdate must be MH or Gibbs")
    }
    if(any(rowSums(y.unk)>1)){#needed if latent histories have >1 capture
      if(!any(rowSums(y.unk>0)>1)){
        probtype="multiple"
      }else{
        probtype="multipletrap"
        
      }
      njs=apply(y.unk,1,function(x){which(x>0)})
      matches=matrix(FALSE,nrow=nUnk,ncol=nUnk)
      diag(matches)=TRUE
      fix=which(duplicated(data.frame(y.unk)))
      
      if(length(fix)>0){
        for(i in 1:length(fix)){
          matches2=which(apply(y.unk,1,function(x){all(x==y.unk[fix[i],])}))
          # matches[fix[i],matches2]=TRUE
          for(j in 1:length(matches2)){
            matches[matches2[j],matches2]=TRUE
          }
        }
      }
      fix=which(duplicated(njs))
      if(length(fix)>0){
        njs=njs[-fix]
      }

      # ll.unk=dbinom(sum(rowSums(y.true)[rowSums(y.true)>0]-1),sum(y.true),theta,log=TRUE)
      
    }else{
      probtype="single"
    }
    if(IDupdate=="MH"){
      # y.true=matrix(0,nrow=M,ncol=J)
      # for(l in 1:nUnk){
      #   nj=which(y.unk[l,]>0)
      #   pii <- lamd[,nj] / sum(lamd[,nj])
      #   ID[l]=sample(1:M,1,prob=pii)
      #   y.true[ID[l],nj]=1
      # }
      S.np <- function(n,m){
        ## Purpose:  Stirling Numbers of the 2-nd kind
        ## 		S^{(m)}_n = number of ways of partitioning a set of
        ##                      $n$ elements into $m$ non-empty subsets
        ## Author: Martin Maechler, Date:  May 28 1992, 23:42
        ## ----------------------------------------------------------------
        ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
        ## Closed Form : p.824 "C."
        ## ----------------------------------------------------------------
        if (0 > m || m > n) stop("'m' must be in 0..n !")
        k <- 0:m
        sig <- rep(c(1,-1)*(-1)^m, length= m+1)
        ga <- gamma(k+1)
        round(sum( sig * k^n /(ga * rev(ga))))
      }
      
      ll.unk=ll.split=ll.order=ll.hyper=rep(0,M)
      idx2=which(rowSums(y.true)>0)
      for(i in idx2){
        maxbreaks=max(sum(y.true[i,])-1,0)
        obs=y.unk[ID==i,]
        
        if(is.matrix(obs)){
          if(nrow(obs)>0){
            obsbreaks=nrow(obs)-1
          }else{
            obsbreaks=0
          }
        }else{
          obsbreaks=0
        }
        ll.unk[i]=dbinom(obsbreaks,maxbreaks,theta,log=TRUE)
        y.sum=sum(y.true[i,])
        if(obsbreaks>0){
          unk.sum=rowSums(y.unk[ID==i,])
          ll.split[i]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks+1),obsbreaks+1),log=TRUE)
        }else{
          unk.sum=sum(y.unk[ID==i,])
          ll.split[i]=0
        }
        if(obsbreaks==0){
          ll.order[i]=0
        }else{
          nj=which(y.true[i,]>0)
          observed <- y.unk[ID==i,nj]
          if(is.null(nrow(observed))){
            observed=matrix(observed,nrow=sum(observed>0))
          }
          # print(unk.sum)
          ll.order[i]=log((factorial(y.sum)/prod(factorial(unk.sum)))/S.np(y.sum,obsbreaks+1)/factorial(obsbreaks+1))
          # ll.order[i]=log((factorial(y.sum)/(prod(factorial(unk.sum))*prod(factorial(table(unk.sum)))))/S.np(y.sum,obsbreaks+1)/factorial(obsbreaks+1))
          ll.hyper[i]=0
          deal=y.true[i,nj]
          if(length(nj)>1){
            for(i2 in 1:(obsbreaks)){
              ll.hyper[i]=ll.hyper[i]+log(prod(choose(deal,observed[i2,]))/choose(sum(deal),unk.sum[i2]))
              deal=deal-observed[i2,]
            }
          }else{
            # ll.hyper[i]=log((factorial(deal)/prod(factorial(unk.sum)))/S.np(deal,obsbreaks+1)/factorial(obsbreaks+1))
            # ll.hyper[i]=log(1/(factorial(deal)/prod(factorial(unk.sum))))
            # ll.hyper[i]=log(1/S.np(deal,obsbreaks+1)/factorial(obsbreaks+1))
            # ll.hyper[i]=log(1/S.np(deal,obsbreaks+1))
            
            # ll.hyper[i]=log(1/choose(deal-1,obsbreaks))
            # ll.hyper[i]=log(1/choose(deal,obsbreaks+1))
            # ll.hyper[i]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks+1),obsbreaks+1),log=TRUE)
            # if(deal>1&(obsbreaks+1)>1){
            #   ll.hyper[i]=log(1/sum(apply(restrictedparts(deal,obsbreaks+1),2,function(x){all(x>0)})))
            # }else{
            #   ll.hyper[i]=0
            # }
            ll.hyper[i]=0
          }
        }
      }
      ll.unk.cand=ll.unk
      ll.split.cand=ll.split
      ll.order.cand=ll.order
      ll.hyper.cand=ll.hyper
    }
    if(all(rowSums(y.unk)==1)){
      theta=1
      uptheta=FALSE
    }else{
      uptheta=TRUE
    }
    
    
    ########################################    
    
    for(iter in 1:niter){
      #Update lam0
      if(obstype=="bernoulli"){
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
          }
        }
      }else{#poisson
        #Update lam0
        llysum=sum(ll.y)
        lam0.cand<- rnorm(1,lam0,proppars$lam0)
        if(lam0.cand > 0){
          lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
          ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp(llycandsum-llysum)){
            lam0<- lam0.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
            llysum=llycandsum
          }
        }
        # Update sigma
        sigma.cand<- rnorm(1,sigma,proppars$sigma)
        if(sigma.cand > 0){
          lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
          ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
          llycandsum=sum(ll.y.cand)
          if(runif(1) < exp((llycandsum)-(llysum))){
            sigma<- sigma.cand
            lamd=lamd.cand
            ll.y=ll.y.cand
          }
        }
      
    }
      
      #Update ID
      if(nUnk>0){
        if(IDupdate=="Gibbs"){
          if(Gibbs=="RC"){#No connections or exclusions
            y.true=y2D
            for(j in trapcaps){
              switchprobs=z*K*lamd[,j]
              switchprobs=switchprobs/sum(switchprobs)
              y.true[,j]=y.true[,j]+rmultinom((n+1):M,njGibbs[j],switchprobs)
            }
            if(obstype=="bernoulli"){
              ll.y= dbinom(y.true,z*pd,K,log=TRUE)
            }else{
              ll.y= dpois(y.true,z*K*lamd,log=TRUE)
            }
            ID=which(rowSums(y.true>0)>0)
            known.vector[(n+1):M]=1*(rowSums(y.true[(n+1):M,])>0)
          }else if(Gibbs=="simple"){#No ID connections across traps
            for(l in 1:nUnk){
              nj=which(y.unk[l,]>0)
              switchprobs=z[(n+1):M]*K*lamd[(n+1):M,nj]
              match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
              if(length(match)!=nUnk){
                match=match[-which(ID[match]==ID[l])]#Can always match with your current cluster members
                legal=rep(TRUE,length(match))
                for(i in 1:length(match)){
                  newcluster=which(ID==ID[match[i]])#who else is linked to the cells you are compatible with?
                  if(length(newcluster)>0){
                    for(j in 1:length(newcluster)){
                      if(constraints[l,newcluster[j]]==0){#are you compatible with their cell, too?
                        legal[i]=FALSE
                      }
                    }
                  }
                }
                rem=match[!legal]
                match=c(l,match[legal])#take out illegal moves and add back in current spot
                if(length(rem)>0){
                  switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
                }
              }
              # if(y.unk[l,nj]>1){
              #   switchprobs=switchprobs^y.unk[l,nj]
              # }
              switchprobs=switchprobs/sum(switchprobs)
              newID=sample((n+1):M,1,prob=switchprobs)

              if(ID[l]!=newID){
                swapped=c(ID[l],newID)
                #update y.true
                y.true[ID[l],]=y.true[ID[l],]-y.unk[l,]
                y.true[newID,]=y.true[newID,]+y.unk[l,]
                ID[l]=newID
                if(obstype=="bernoulli"){
                  ll.y[swapped,]= dbinom(y.true[swapped,],K,pd[swapped,],log=TRUE)
                }else{
                  ll.y[swapped,]= dpois(y.true[swapped,],K*lamd[swapped,],log=TRUE)
                }
              }
            }
            known.vector[(n+1):M]=1*(rowSums(y.true[(n+1):M,])>0)
          }

        }else{ # MH update

          if(1==2){
          up=sample(1:nUnk,nswap,replace=FALSE)
          for(l in up){
            #traps where this sample caught
            nj=which(y.unk[l,]>0)
            # possible=(n+1):M
            # possible=possible[z==1]
            # if(length(possible)>0){#Check to see if constraints match
            #   thiscluster=which(ID==ID[l])
            #   if(length(thiscluster)>1){
            #     thiscluster=thiscluster[-which(rowSums(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
            #   }else{
            #     thiscluster=thiscluster[-which(sum(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
            #   }
            #   if(length(thiscluster)>0){
            #     legal=rep(TRUE,length(possible))
            #     #Check to see if you can swap into these clusters
            #     for(i in 1:length(possible)){
            #       check=which(ID==possible[i])#Who else is currently assigned this possible new ID?
            #       if(length(check)>0){#if false, no ID assigned to this guy and legal stays true
            #         for(k in 1:length(thiscluster)){#loop through guys in focal cluster
            #           for(j in 1:length(check)){#check against guys in possible cluster
            #             if(constraints[thiscluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
            #               legal[i]=FALSE
            #             }
            #           }
            #         }
            #       }
            #     }
            #     possible=possible[legal]
            #   }
            #   if(length(possible)>0){#Can update
                newID=ID
                if(probtype=="single"|probtype=="multiple"){
                  njprobs=lamd[,nj]*z #z=1 for possible guys
                  njprobs=njprobs/sum(njprobs)
                  newID[l]=sample(1:M,1,prob=njprobs)
                }else{
                  # njprobs=matrix(NA,nrow=M,ncol=length(nj))
                  # for(j in 1:length(nj)){
                  #   njprobs[,j]=lamd[,nj[j]]*z
                  #   njprobs[,j]=njprobs[,j]/sum(njprobs[,j])
                  # }
                  # njprobs=apply(njprobs,1,prod)
                  dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
                  possible<- which(dv < swap.tol)+n
                  possible=possible[z[possible]==1]
                  possible=possible[-which(possible==ID[l])]
                  if(length(possible)==0)next
                  jump.probability<- 1/length(possible)  #
                  if(length(possible)>1){
                    newID[l]=sample(possible,1)
                  }else{
                    newID[l]=possible
                  }
                }
                
                # njprobs=njprobs^y.unk[l,nj]
                # njprobs=njprobs/sum(njprobs)
                # if(length(possible)>1){
                  # newID[l]=sample((n+1):M,1,prob=njprobs)
                # }
               
                # choose1=which(possible==newID[l])
                # trash=(n+1):M
                # trash=trash[z==1]
                # backcluster=which(newID==newID[l])#which guys are in the new ID cluster?
                # if(length(backcluster)>1){
                #   backcluster=backcluster[-which(rowSums(constraints[backcluster,])!=nUnk)]#Don't need to check if no constraints for this sc
                # }else{
                #   backcluster=backcluster[-which(sum(constraints[backcluster,])==nUnk)]
                # }
                # # Check to see if you can swap into these clusters
                # if(length(backcluster)>0){
                #   legal=rep(TRUE,length(trash))
                #   for(i in 1:length(trash)){
                #     check=which(newID==trash[i])
                #     for(k in 1:length(backcluster)){#loop through guys in back cluster
                #       if(length(check)>0){#if false, no ID assigned to this back guy and legal stays true
                #         for(j in 1:length(check)){#check against back guys in proposed cluster
                #           if(constraints[backcluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                #             legal[i]=FALSE
                #           }
                #         }
                #       }
                #     }
                #   }
                #   trash=trash[legal]
                # }
                swapped=c(ID[l],newID[l])#order swap.out then swap.in
                cluster1guys=which(newID==swapped[1])
                cluster2guys=which(newID==swapped[2])
                y.cand=array(0,dim=c(2,J))
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
                # if(all(possible!=trash)){#are back probs going to be different from forward?
                #   if(length(nj)==1){
                #     njbackprobs=lamd[trash,nj] #z=1 for possible guys
                #   }else{
                #     njbackprobs=matrix(NA,nrow=length(trash),ncol=length(nj))
                #     for(j in 1:length(nj)){
                #       njbackprobs[,j]=lamd[trash,nj[j]]
                #       njbackprobs[,j]=njbackprobs[,j]/sum(njbackprobs[,j])
                #     }
                #     njbackprobs=apply(njbackprobs,1,prod)
                #   }
                #   njbackprobs=njbackprobs/sum(njbackprobs)
                # }else{
                #   njbackprobs=njprobs
                # }
                # choose2=which(trash==ID[l])
                # focalprob=1/nUnk
                # njbackprobs=njprobs
                # focalprob=1/sum(y.unk[,nj]>0&ID==ID[l])
                y.cand2=y.true
                y.cand2[swapped,]=y.cand
                if(probtype=="single"){
                  focalprob=1/nUnk
                  focalbackprob=1/sum(y.unk[,nj]>0&newID==newID[l])
                }else if(probtype=="multiple"){
                  focalbackprob=sum(y.unk[,nj]==y.unk[l,nj]&newID==newID[l])/sum(y.unk[,nj]>0&newID==newID[l])
                  focalprob=1/nUnk
                }else{
                  # focalbackprob=1/sum(newID==newID[l])
                  if(length(nj)>1){
                    #forward probs
                    # b=apply(y.unk[,nj],1,function(x){all(x>0)})&rowSums(y.unk[,-nj])==0#samples at all of these traps and not others
                    # focalprob=1/nUnk
                    # focalprob=1/sum(b&ID==ID[l])
                    a=apply(y.unk[,nj],1,function(x){all(x==y.unk[l,nj])})
                    focalprob=1/sum(a&rowSums(y.unk[,-nj])==0&ID==ID[l])#exact matches
                  
                    
                    # focalbackprob=sum(a&rowSums(y.unk[,-nj])==0&newID==newID[l])/sum(apply(y.unk[,nj],1,function(x){all(x>0)})&rowSums(y.unk[,-nj])==0&newID==newID[l])
                    # focalbackprob=sum(a&newID==newID[l])/sum(apply(y.unk[,nj],1,function(x){all(x>0)})&newID==newID[l])
                    focalbackprob=1/sum(a&rowSums(y.unk[,-nj])==0&newID==newID[l])#exact matches
                    # focalbackprob=1/sum(a&newID==newID[l])
                 
                  }else{
                    #forward probs
                    # b=y.unk[,nj]>0&rowSums(y.unk[,-nj])==0#samples at this trap and not others
                    # focalprob=sum(b)/nUnk*1/length(unique(ID[b]))*1/sum(ID[b]==ID[l])
                    #backwards probs
                    # focalprob=1/sum(b&ID==ID[l])
                    focalprob=1/sum(y.unk[,nj]==y.unk[l,nj]&rowSums(y.unk[,-nj])==0&ID==ID[l])#exact matches
                    focalbackprob=1/sum(y.unk[,nj]==y.unk[l,nj]&rowSums(y.unk[,-nj])==0&newID==newID[l])#exact matches
                    # focalbackprob=1/sum(y.unk[,nj]==y.unk[l,nj]&newID==newID[l])
                     
                  }
                }
                
                
                ll.y.cand[swapped,]=dpois(y.cand,K*lamd[swapped,],log=TRUE)
                y.cand2=y.true
                y.cand2[swapped,]=y.cand
                maxbreaks=rowSums(y.cand2[swapped,])-1
                maxbreaks[maxbreaks<0]=0
                obsbreaks=rep(NA,2)
                obs1=y.unk[newID==swapped[1],]
                obs2=y.unk[newID==swapped[2],]
                if(is.matrix(obs1)){
                  if(nrow(obs1)>0){
                    obsbreaks[1]=max(0,nrow(obs1)-1)
                  }else{
                    obsbreaks[1]=0
                  }
                }else{
                  obsbreaks[1]=0
                }
                if(is.matrix(obs2)){
                  if(nrow(obs2)>0){
                    obsbreaks[2]=max(0,nrow(obs2)-1)
                  }else{
                    obsbreaks[2]=0
                  }
                  
                }else{
                  obsbreaks[2]=0
                }
                
                ll.unk.cand[swapped]=dbinom(obsbreaks,maxbreaks,theta,log=TRUE)
                dv<-  sqrt( (s[newID[l],1]- s[(n+1):M,1])^2 + (s[newID[l],2] - s[(n+1):M,2])^2 )
                trash<-   which(dv < swap.tol)+n
                trash=trash[z[trash]==1]
                trash=trash[-which(trash==newID[l])]
                back.probability=1/length(trash)
                
                y.sum=sum(y.cand2[ID[l],])
                if(y.sum==0){
                  ll.split.cand[ID[l]]=0
                  ll.order.cand[ID[l]]=0
                }else{
                  if(obsbreaks[1]>0){
                    unk.sum=rowSums(y.unk[newID==ID[l],])
                    ll.split.cand[ID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[1]+1),obsbreaks[1]+1),log=TRUE)
                    ll.order.cand[ID[l]]=log((factorial(y.sum)/prod(factorial(unk.sum)))/S.np(y.sum,obsbreaks[1]+1)/factorial(obsbreaks[1]+1))
                  }else{
                    ll.split.cand[ID[l]]=0
                    ll.order.cand[ID[l]]=0
                  }
                }
                #update new guy
                y.sum=sum(y.cand2[newID[l],])
                if(y.sum==0){#Should never happen?
                  ll.split.cand[newID[l]]=0
                  ll.order.cand[newID[l]]=0
                }else{
                  if(obsbreaks[2]>0){
                    unk.sum=rowSums(y.unk[newID==newID[l],])
                    ll.split.cand[newID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[2]+1),obsbreaks[2]+1),log=TRUE)
                    ll.order.cand[newID[l]]=log((factorial(y.sum)/prod(factorial(unk.sum)))/S.np(y.sum,obsbreaks[2]+1)/factorial(obsbreaks[2]+1))
               
                  }else{
                    ll.split.cand[newID[l]]=0
                    ll.order.cand[newID[l]]=0
                  }
                }

               
                ll.y.cand[swapped,]=dpois(y.cand,K*lamd[swapped,],log=TRUE)
                # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.split.cand[swapped])+sum(ll.order.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.split[swapped])+sum(ll.order[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.order.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.order[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  
                # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.split.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.split[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  y.true[swapped,]=y.cand
                  ll.y[swapped,]=ll.y.cand[swapped,]
                  ll.unk[swapped]=ll.unk.cand[swapped]
                  ll.split[swapped]=ll.split.cand[swapped]
                  ll.order[swapped]=ll.order.cand[swapped]
                  ID[l]=newID[l]
                  known.vector[swapped]=1*(rowSums(y.cand)>0)
                }
               
          }
          }else{
            if(probtype!="multipletrap"){
              up=which(colSums(y.unk)>0)
              for(nj in up){
                #pick an ID to update
                ids=which(y.true[,nj]>0) #who is assigned to samples at this trap?
                if(length(ids)>1){
                  IDpick=sample(ids,1) #pick one ID
                }else{
                  IDpick=ids
                }
                samples=which(ID==IDpick&y.unk[,nj]>0) #pick one sample|ID
                if(length(samples)>1){
                  l=sample(samples,1) #which sample do we choose?
                }else{
                  l=samples
                }
                # possible=(n+1):M
                # possible=possible[z==1]
                # if(length(possible)>0){#Check to see if constraints match
                #   thiscluster=which(ID==ID[l])
                #   if(length(thiscluster)>1){
                #     thiscluster=thiscluster[-which(rowSums(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
                #   }else{
                #     thiscluster=thiscluster[-which(sum(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
                #   }
                #   if(length(thiscluster)>0){
                #     legal=rep(TRUE,length(possible))
                #     #Check to see if you can swap into these clusters
                #     for(i in 1:length(possible)){
                #       check=which(ID==possible[i])#Who else is currently assigned this possible new ID?
                #       if(length(check)>0){#if false, no ID assigned to this guy and legal stays true
                #         for(k in 1:length(thiscluster)){#loop through guys in focal cluster
                #           for(j in 1:length(check)){#check against guys in possible cluster
                #             if(constraints[thiscluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                #               legal[i]=FALSE
                #             }
                #           }
                #         }
                #       }
                #     }
                #     possible=possible[legal]
                #   }
                #   if(length(possible)>0){#Can update
                newID=ID
                # if(length(nj)==1){
                #   njprobs=lamd[,nj]*z #z=1 for possible guys
                #   # njprobs=lamd[possible,nj] #z=1 for possible guys
                #   
                # }else{
                #   njprobs=matrix(NA,nrow=length(possible),ncol=length(nj))
                #   for(j in 1:length(nj)){
                #     njprobs[,j]=lamd[possible,nj[j]]
                #     njprobs[,j]=njprobs[,j]/sum(njprobs[,j])
                #   }
                #   njprobs=apply(njprobs,1,prod)
                # }
                # njprobs=njprobs/sum(njprobs)
                newID[l]=sample.int(length(njprobs),1,prob=njprobs)
                dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
                possible<- which(dv < swap.tol)+n
                possible=possible[z[possible]==1]
                jump.probability<- 1/length(possible)  #
                
                if(length(possible)>1){
                  # newID[l]=sample(possible,1,prob=njprobs)
                  newID[l]=sample.int(length(njprobs),1,prob=njprobs)
                }else{
                  newID[l]=possible
                }
                if(ID[l]!=newID[l]){
                  # choose1=which(possible==newID[l])
                  # choose1=newID[l]
                  # trash=(n+1):M
                  # trash=trash[z==1]
                  # backcluster=which(newID==newID[l])#which guys are in the new ID cluster?
                  # if(length(backcluster)>1){
                  #   backcluster=backcluster[-which(rowSums(constraints[backcluster,])!=nUnk)]#Don't need to check if no constraints for this sc
                  # }else{
                  #   backcluster=backcluster[-which(sum(constraints[backcluster,])==nUnk)]
                  # }
                  # # Check to see if you can swap into these clusters
                  # if(length(backcluster)>0){
                  #   legal=rep(TRUE,length(trash))
                  #   for(i in 1:length(trash)){
                  #     check=which(newID==trash[i])
                  #     for(k in 1:length(backcluster)){#loop through guys in back cluster
                  #       if(length(check)>0){#if false, no ID assigned to this back guy and legal stays true
                  #         for(j in 1:length(check)){#check against back guys in proposed cluster
                  #           if(constraints[backcluster[k],check[j]]==0){#if any members of the cluster are inconsistent, illegal move
                  #             legal[i]=FALSE
                  #           }
                  #         }
                  #       }
                  #     }
                  #   }
                  #   trash=trash[legal]
                  # }
                  swapped=c(ID[l],newID[l])#order swap.out then swap.in
                  cluster1guys=which(newID==swapped[1])
                  cluster2guys=which(newID==swapped[2])
                  y.cand=array(0,dim=c(2,J))
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
                  # if(all(possible!=trash)){#are back probs going to be different from forward?
                  #   if(length(nj)==1){
                  #     njbackprobs=lamd[trash,nj] #z=1 for possible guys
                  #   }else{
                  #     njbackprobs=matrix(NA,nrow=length(trash),ncol=length(nj))
                  #     for(j in 1:length(nj)){
                  #       njbackprobs[,j]=lamd[trash,nj[j]]
                  #       njbackprobs[,j]=njbackprobs[,j]/sum(njbackprobs[,j])
                  #     }
                  #     njbackprobs=apply(njbackprobs,1,prod)
                  #   }
                  #   njbackprobs=njbackprobs/sum(njbackprobs)
                  # }else{
                  #   njbackprobs=njprobs
                  # }
                  njbackprobs=njprobs
                  # choose2=which(trash==ID[l])
                  # choose2=ID[l]
                  #backprobs
                  y.cand2=y.true
                  y.cand2[swapped,]=y.cand
                  idsback=which(y.cand2[,nj]>0) #who is assigned to samples at this trap?
                  # IDpickback=newID[l] #must pick this guy to swap back
                  samplesback=which(newID==newID[l]&y.unk[,nj]>0) #pick
                  focalprob=(1/length(ids))*(sum(y.unk[ID==ID[l]&y.unk[,nj]>0,nj]==y.unk[l,nj])/length(samples))#*(1/length(samples))
                  focalbackprob=(1/length(idsback))*(sum(y.unk[newID==newID[l]&y.unk[,nj]>0,nj]==y.unk[l,nj])/length(samplesback))#(1/length(samplesback))
                  # p1=log(focalprob*njprobs[newID[l]])
                  # p2=log(focalbackprob*njbackprobs[ID[l]])
                  maxbreaks=rowSums(y.cand2[swapped,])-1
                  maxbreaks[maxbreaks<0]=0
                  obsbreaks=rep(NA,2)
                  obs1=y.unk[newID==swapped[1],]
                  obs2=y.unk[newID==swapped[2],]
                  if(is.matrix(obs1)){
                    if(nrow(obs1)>0){
                      obsbreaks[1]=max(0,nrow(obs1)-1)
                    }else{
                      obsbreaks[1]=0
                    }
                  }else{
                    obsbreaks[1]=0
                  }
                  if(is.matrix(obs2)){
                    if(nrow(obs2)>0){
                      obsbreaks[2]=max(0,nrow(obs2)-1)
                    }else{
                      obsbreaks[2]=0
                    }
                    
                  }else{
                    obsbreaks[2]=0
                  }
                  #
                  ll.unk.cand[swapped]=dbinom(obsbreaks,maxbreaks,theta,log=TRUE)

                  
                  ll.y.cand[swapped,]=dpois(y.cand,K*lamd[swapped,],log=TRUE)
                    if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped]))-(sum(ll.y[swapped,])+sum(ll.unk[swapped])))*((njbackprobs[ID[l]]*focalbackprob)/(njprobs[newID[l]]*focalprob))){
                    y.true[swapped,]=y.cand
                    ll.y[swapped,]=ll.y.cand[swapped,]
                    ll.unk[swapped]=ll.unk.cand[swapped]
                    ID[l]=newID[l]
                    known.vector[swapped]=1*(rowSums(y.cand)>0)
                    ll.unk[swapped]=ll.unk.cand[swapped]
                  }
                }
              }
            }else{
              for(j in 1:length(njs)){
                nj=njs[j][[1]]
                #pick an ID to update
                if(length(nj)==1){
                  ids=which(y.true[,nj]>0) #who is assigned to samples at this trap?
                  #remove ids that don't have any samples we can move
                  rem=c()
                  for(i in 1:length(ids)){
                    if(!any(ID==ids[i]&rowSums(y.unk[,-nj]>0)==0&y.unk[,nj]>0)){
                      rem=c(rem,i)
                    }
                  }
                  if(length(rem)>0){
                    ids=ids[-rem]
                  }
                }else{
                  ids=which(rowSums(y.true[,nj]>0)==length(nj))
                  #remove ids involving other traps
                  rem=c()
                  for(i in 1:length(ids)){
                    if(!any(ID==ids[i]&rowSums(y.unk[,nj]>0)==length(nj)&rowSums(y.unk[,-nj]>0)==0)){
                      rem=c(rem,i)
                    }
                  }
                  if(length(rem)>0){
                    ids=ids[-rem]
                  }
                }
                if(length(ids)>1){
                  IDpick=sample(ids,1) #pick one ID
                }else{
                  IDpick=ids
                }
                if(length(nj)==1){
                  samples=which(ID==IDpick&y.unk[,nj]>0&rowSums(y.unk[,-nj]>0)==0) #pick one sample|ID
                }else{
                  samples=which(ID==IDpick&rowSums(y.unk[,nj]>0)==length(nj)&rowSums(y.unk[,-nj]>0)==0) #pick one sample|ID
                }
                if(length(samples)>1){
                  l=sample(samples,1) #which sample do we choose?
                }else{
                  l=samples
                }
                # if(length(nj)==1){
                #   njprobs=lamd[,nj]*z #z=1 for possible guys
                #   # njprobs=lamd[possible,nj] #z=1 for possible guys
                #   
                # }else{
                #   njprobs=matrix(NA,nrow=M,ncol=length(nj))
                #   for(j in 1:length(nj)){
                #     njprobs[,j]=z*lamd[,nj[j]]
                #     njprobs[,j]=njprobs[,j]/sum(njprobs[,j])
                #   }
                #   njprobs=apply(njprobs,1,prod)
                # }
                # njprobs=njprobs/sum(njprobs)
                # 
                # newID[l]=sample.int(length(njprobs),1,prob=njprobs)
                
                dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
                possible<- which(dv < swap.tol)+n
                possible=possible[z[possible]==1]
                if(length(possible)>0){
                  possible=possible[-which(possible==ID[l])]
                }
                if(length(possible)>0){
                jump.probability<- 1/length(possible)  
                newID=ID
                if(length(possible)>1){
                  newID[l]=sample(possible,1)
                }else{
                  newID[l]=possible
                }

                  swapped=c(ID[l],newID[l])#order swap.out then swap.in
                  cluster1guys=which(newID==swapped[1])
                  cluster2guys=which(newID==swapped[2])
                  y.cand=array(0,dim=c(2,J))
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
                  # njbackprobs=njprobs
                  y.cand2=y.true
                  y.cand2[swapped,]=y.cand
                  if(length(nj)==1){
                    idsback=which(y.cand2[,nj]>0) #who is assigned to samples at this trap?
                    #remove ids that don't have at least one sample we can move
                    rem=c()
                    for(i in 1:length(idsback)){
                      # if(!any(newID==idsback[i]&y.unk[,nj]>0&rowSums(y.unk[,-nj]>0)==0)){
                      if(!any(newID==idsback[i]&matches[l,])){#exact match
                        rem=c(rem,i)
                      }
                    }
                    if(length(rem)>0){
                      idsback=idsback[-rem]
                    }
                    # samplesback=which(newID==newID[l]&y.unk[,nj]>0&rowSums(y.unk[,-nj])==0) #pick
                    samplesback=which(newID==newID[l]&matches[l,]) #exact matches
                    focalprob=(1/length(ids))*(1/length(samples))
                    focalbackprob=(1/length(idsback))*(1/length(samplesback))
                  }else{
                    idsback=which(rowSums(y.cand2[,nj]>0)==length(nj)) #who is assigned to samples at this trap?
                    #remove idsback involving other traps
                    rem=c()
                    for(i in 1:length(idsback)){
                      # if(!any(newID==idsback[i]&rowSums(y.unk[,nj]>0)==length(nj)&rowSums(y.unk[,-nj]>0)==0)){
                      if(!any(newID==idsback[i]&matches[l,])){#exact matches
                        rem=c(rem,i)
                      }
                    }
                    if(length(rem)>0){
                      idsback=idsback[-rem]
                    }
                    # samplesback=which(newID==newID[l]&rowSums(y.unk[,nj]>0)==length(nj)&rowSums(y.unk[,-nj]>0)==0)
                    samplesback=which(newID==newID[l]&matches[l,]) #exact matches
                    focalprob=(1/length(ids))*(1/length(samples))
                    focalbackprob=(1/length(idsback))*(1/length(samplesback))
                  }
                  maxbreaks=rowSums(y.cand2[swapped,])-1
                  maxbreaks[maxbreaks<0]=0
                  obsbreaks=rep(NA,2)
                  obs1=y.unk[which(newID==swapped[1]),]
                  obs2=y.unk[which(newID==swapped[2]),]
                  if(is.matrix(obs1)){
                    if(nrow(obs1)>0){
                      obsbreaks[1]=max(0,nrow(obs1)-1)
                    }else{
                      obsbreaks[1]=0
                    }
                  }else{
                    obsbreaks[1]=0
                  }
                  if(is.matrix(obs2)){
                    if(nrow(obs2)>0){
                      obsbreaks[2]=max(0,nrow(obs2)-1)
                    }else{
                      obsbreaks[2]=0
                    }
                    
                  }else{
                    obsbreaks[2]=0
                  }
                  
                  ll.unk.cand[swapped]=dbinom(obsbreaks,maxbreaks,theta,log=TRUE)
                  dv<-  sqrt( (s[newID[l],1]- s[(n+1):M,1])^2 + (s[newID[l],2] - s[(n+1):M,2])^2 )
                  trash<-   which(dv < swap.tol)+n
                  trash=trash[z[trash]==1]
                  trash=trash[-which(trash==newID[l])]
                  back.probability=1/length(trash)

                  #update old guy
                  y.sum=sum(y.cand2[ID[l],])
                  if(y.sum==0){
                    ll.split.cand[ID[l]]=0
                    ll.order.cand[ID[l]]=0
                    ll.hyper.cand[ID[l]]=0
                  }else{
                    if(obsbreaks[1]>0){
                      unk.sum=rowSums(y.unk[newID==ID[l],])
                      ll.split.cand[ID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[1]+1),obsbreaks[1]+1),log=TRUE)
                      ll.order.cand[ID[l]]=log((factorial(y.sum)/prod(factorial(unk.sum)))/S.np(y.sum,obsbreaks[1]+1)/factorial(obsbreaks[1]+1))
                      # ll.order.cand[ID[l]]=log((factorial(y.sum)/(prod(factorial(unk.sum))* prod(factorial(table(unk.sum)))))/S.np(y.sum,obsbreaks[1]+1)/factorial(obsbreaks[1]+1))
                     
                      nj=which(y.cand2[ID[l],]>0)
                      observed=obs1[,nj]
                      if(!is.matrix(observed)){
                        observed=matrix(observed,nrow=sum(observed>0))
                      }
                      deal=y.cand2[ID[l],nj]
                      ll.hyper.cand[ID[l]]=0
                      if(length(nj)>1){
                        for(i2 in 1:(obsbreaks[1])){
                          ll.hyper.cand[ID[l]]=ll.hyper.cand[ID[l]]+log(prod(choose(deal,observed[i2,]))/choose(sum(deal),unk.sum[i2]))
                          deal=deal-observed[i2,]
                        }
                      }else{
                        ll.hyper.cand[ID[l]]=log((factorial(deal)/prod(factorial(unk.sum)))/S.np(deal,obsbreaks[1]+1)/factorial(obsbreaks[1]+1))
                        # ll.hyper[ID[l]]=log(1/(factorial(deal)/prod(factorial(unk.sum))))
                        # ll.hyper[ID[l]]=log(1/S.np(deal,obsbreaks[1]+1)/factorial(obsbreaks[1]+1))
                        # ll.hyper[ID[l]]=log(1/S.np(deal,obsbreaks[1]+1))
                        
                        # ll.hyper[ID[l]]=log(1/choose(deal-1,obsbreaks[1]))
                        # ll.hyper[ID[l]]=log(1/choose(deal,obsbreaks[1]+1))
                        # ll.hyper[ID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[1]+1),obsbreaks[1]+1),log=TRUE)
                        
                        # if(deal>1&(obsbreaks[1]+1)>1){
                        #   ll.hyper[ID[l]]=log(1/sum(apply(restrictedparts(deal,obsbreaks[1]+1),2,function(x){all(x>0)})))
                        # }else{
                        #   ll.hyper[ID[l]]=0
                        # }
                        ll.hyper.cand[ID[l]]=0
                        
                      }
                    }else{
                      ll.split.cand[ID[l]]=0
                      ll.order.cand[ID[l]]=0
                      ll.hyper.cand[ID[l]]=0
                    }
                  }
              
                  
                  #update new guy
                  y.sum=sum(y.cand2[newID[l],])
                  if(y.sum==0){#Should never happen?
                    ll.split.cand[newID[l]]=0
                    ll.order.cand[newID[l]]=0
                    ll.hyper.cand[newID[l]]=0
                  }else{
                    if(obsbreaks[2]>0){
                      unk.sum=rowSums(y.unk[newID==newID[l],])
                      ll.split.cand[newID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[2]+1),obsbreaks[2]+1),log=TRUE)
                      ll.order.cand[newID[l]]=log((factorial(y.sum)/prod(factorial(unk.sum)))/S.np(y.sum,obsbreaks[2]+1)/factorial(obsbreaks[2]+1))
                      # ll.order.cand[newID[l]]=log((factorial(y.sum)/(prod(factorial(unk.sum))*prod(factorial(table(unk.sum)))))/S.np(y.sum,obsbreaks[2]+1)/factorial(obsbreaks[2]+1))
                      
                      nj=which(y.cand2[newID[l],]>0)
                      observed=obs2[,nj]
                      if(!is.matrix(observed)){
                        observed=matrix(observed,nrow=sum(observed>0))
                      }
                      deal=y.cand2[newID[l],nj]
                      ll.hyper.cand[newID[l]]=0
                      if(length(nj)>1){
                        for(i2 in 1:(obsbreaks[2])){
                          ll.hyper.cand[newID[l]]=ll.hyper.cand[newID[l]]+log(prod(choose(deal,observed[i2,]))/choose(sum(deal),unk.sum[i2]))
                          deal=deal-observed[i2,]
                        }
                      }else{
                        # ll.hyper.cand[newID[l]]=log((factorial(deal)/prod(factorial(unk.sum)))/S.np(deal,obsbreaks[2]+1)/factorial(obsbreaks[2]+1))
                        # ll.hyper[newID[l]]=log(1/(factorial(deal)/prod(factorial(unk.sum))))
                        # ll.hyper[newID[l]]=log(1/S.np(deal,obsbreaks[2]+1)/factorial(obsbreaks[2]+1))
                        # ll.hyper[newID[l]]=log(1/S.np(deal,obsbreaks[2]+1))
                        # ll.hyper[newID[l]]=log(1/choose(deal-1,obsbreaks[2]))
                        
                        # ll.hyper[newID[l]]=log(1/choose(deal,obsbreaks[2]+1))
                        # ll.hyper[newID[l]]=dmultinom(unk.sum-1,y.sum-length(unk.sum),prob=rep(1/(obsbreaks[2]+1),obsbreaks[2]+1),log=TRUE)
                        # if(deal>1&(obsbreaks[2]+1)>1){
                        #   ll.hyper[newID[l]]=log(1/sum(apply(restrictedparts(deal,obsbreaks[2]+1),2,function(x){all(x>0)})))
                        # }else{
                        #   ll.hyper[newID[l]]=0
                        # }
                        ll.hyper[newID[l]]=0
                      }
                    }else{
                      ll.split.cand[newID[l]]=0
                      ll.order.cand[newID[l]]=0
                      ll.hyper.cand[newID[l]]=0
                    }
                  }
                  ll.y.cand[swapped,]=dpois(y.cand,K*lamd[swapped,],log=TRUE)
                  if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.hyper.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.hyper[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.hyper.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.hyper[swapped]))))*(back.probability/jump.probability)){
                    # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.split.cand[swapped])+sum(ll.hyper.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.split[swapped])+sum(ll.hyper[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                      
                  
                  
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.order.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.order[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                    
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.order.cand[swapped])+sum(ll.hyper.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.order[swapped])+sum(ll.hyper[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.split.cand[swapped])+sum(ll.order.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.split[swapped])+sum(ll.order[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.unk.cand[swapped])+sum(ll.split.cand[swapped]))-((sum(ll.y[swapped,])+sum(ll.unk[swapped])+sum(ll.split[swapped]))))*(back.probability/jump.probability)*(focalbackprob/focalprob)){
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,]))-(sum(ll.y[swapped,])))*(njbackprobs[ID[l]]/njprobs[newID[l]])*(focalbackprob/focalprob)){
                    y.true[swapped,]=y.cand
                    ll.y[swapped,]=ll.y.cand[swapped,]
                    ll.unk[swapped]=ll.unk.cand[swapped]
                    ll.split[swapped]=ll.split.cand[swapped]
                    ll.order[swapped]=ll.order.cand[swapped]
                    ll.hyper[swapped]=ll.hyper.cand[swapped]
                    ID[l]=newID[l]
                    known.vector[swapped]=1*(rowSums(y.cand)>0)
                  }
              }
              
            }
          }
          }
        }
      }
      if(uptheta){
        
        theta.cand <- rnorm(1, theta, proppars$theta)
        if(theta.cand >= 0 & theta.cand <= 1) {
          ll.unk.cand=ll.unk*0
          idx2=which(rowSums(y.true)>0)
          for(i in idx2){
            maxbreaks=max(sum(y.true[i,])-1,0)
            obs=y.unk[ID==i,]

            if(is.matrix(obs)){
              if(nrow(obs)>0){
              obsbreaks=nrow(obs)-1
              }else{
                obsbreaks=0
              }
            }else{
              obsbreaks=0
            }
            ll.unk.cand[i]=dbinom(obsbreaks,maxbreaks,theta.cand,log=TRUE)

          }
          # for(l in 1:nUnk){
          #   ll.unk.cand[l]=sum(dbinom(sum(y.unk[l,]),sum(y.true[ID[l],]),theta.cand,log=TRUE))
          # }
          if(runif(1) < exp(sum(ll.unk.cand) -  sum(ll.unk))) {
            ll.unk=ll.unk.cand
            theta = theta.cand
          }
        }
      }
      # if(!all(is.finite(ll.unk.cand)))stop("!")
    
    
      ## probability of not being captured in a trap AT ALL
      if(obstype=="poisson"){
        pd=1-exp(-lamd)
      }
      pbar=(1-pd)^K
      prob0<- exp(rowSums(log(pbar)))
      fc<- prob0*psi/(prob0*psi + 1-psi)
      z[known.vector==0]<- rbinom(sum(known.vector ==0), 1, fc[known.vector==0])
      if(obstype=="bernoulli"){
        ll.y= dbinom(y.true,K,pd*z,log=TRUE)
      }else{
        ll.y= dpois(y.true,K*lamd*z,log=TRUE)
      }

      # zUps <- 0
      # seen <- rowSums(y) > 0 #(y>0, 1, any)
      # for(i in 1:M) {
      #   if(seen[i])
      #     next
      #   zcand <- z
      #   zcand[i] <- ifelse(z[i]==0, 1, 0)
      #     ll <- sum(dpois(y.true[i,], K*lamd[i,]*z[i], log=TRUE))
      #     llcand <- sum(dpois(y.true[i,], K*lamd[i,]*zcand[i], log=TRUE))
      # 
      #   prior <- dbinom(z[i], 1, psi, log=TRUE)
      #   prior.cand <- dbinom(zcand[i], 1, psi, log=TRUE)
      #   if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
      #     z[i] <- zcand[i]
      #     zUps <- zUps+1
      #   }
      # }
      
      psi=rbeta(1,1+sum(z),1+M-sum(z))
      
      
      
      
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
          if(obstype=="bernoulli"){
            pd.cand[i,]=1-exp(-lamd.cand[i,])
            ll.y.cand[i,]= dbinom(y.true[i,],K,pd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              pd[i,]=pd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }else{#poisson
            ll.y.cand[i,]= dpois(y.true[i,],K*lamd.cand[i,]*z[i],log=TRUE)
            if (runif(1) < exp(sum(ll.y.cand[i,]) - sum(ll.y[i,]))) {
              s[i,]=Scand
              D[i,]=dtmp
              lamd[i,]=lamd.cand[i,]
              ll.y[i,]=ll.y.cand[i,]
            }
          }
        }
      }
      #Do we record output on this iteration?
      if(iter>nburn&iter%%nthin==0){
        sxout[idx,]<- s[,1]
        syout[idx,]<- s[,2]
        zout[idx,]<- z
        if(IDupdate=="Gibbs"){
          if(Gibbs!="RC"){
            IDout[idx,]=ID
          }
        }else{
          IDout[idx,]=ID
        }
        out[idx,]<- c(lam0,sigma ,sum(z),length(unique(ID)),theta)
        # print(out[idx,])
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else{ 
      list(out=out)
    }
  }
