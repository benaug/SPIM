genSPIMmcmc2D <-
  function(data,niter=2400,nburn=1200, nthin=5, M = 200,K=NA, inits=inits,obstype="bernoulli",
           proppars=list(lam0=0.05,sigma=0.1,sx=0.2,sy=0.2),keepACs=TRUE,swap.tol=1,swap=10,IDupdate="Gibbs"){
    library(abind)
    y<-data$yID
    y.unk<-data$yUnkobs
    X<-as.matrix(data$X)
    J<-nrow(X)
    n<- dim(y)[1]
    K=data$K
    if(nrow(y.unk)<swap){
      swap=nrow(y.unk)
    }
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
    out<-matrix(NA,nrow=nstore,ncol=4)
    dimnames(out)<-list(NULL,c("lam0","sigma","N","clusters"))
    sxout<- syout<- zout<-matrix(NA,nrow=nstore,ncol=M)
    IDout=matrix(NA,nrow=nstore,ncol=length(ID))
    idx=1 #for storing output not recorded every iteration
    
    #collapse data to 2D
    y.true=apply(y.true,c(1,2),sum)
    y.unk=apply(y.unk,c(1,2),sum)
    
    
    D<- e2dist(s, X)
    lamd=lamd.cand=lam0*exp(-D*D/(2*sigma*sigma))
    if(obstype=="bernoulli"){
      pd=pd.cand=1-exp(-lamd)
      ll.y=ll.y.cand= dbinom(y.true,K,pd*z,log=TRUE)
    }else{
      ll.y=ll.y.cand= dpois(y.true,K*lamd*z,log=TRUE)
    }
    if(nUnk>0&IDupdate=="Gibbs"){
      ntraps=rowSums(y.unk>0)
      if(any(ntraps>1)){#at least one connection across traps
        Gibbs="complex"
        trapcaps=apply(y.unk,1,function(x){which(x>0)})
        overlap=vector("list",nUnk)
        trapcaps2=vector("list",nUnk)
        for(l in 1:nUnk){
          nj=which(y.unk[l,]>0)
          if(length(nj)>0){
            #guys you have overlaps with, but not complete
            overlap[[l]]=which(unlist(lapply(trapcaps,function(x){any(x%in%nj)&!all(x%in%nj)})))
            if(length(overlap[[l]])>0){
              trapcaps2[[l]]=vector("list",length(overlap[[l]]))
              for(i in 1:length(overlap[[l]])){
                trapcaps2[[l]][[i]]=which(y.unk[overlap[[l]][i],]>0)
                trapcaps2[[l]][[i]]=trapcaps2[[l]][[i]][-which(trapcaps2[[l]][[i]]%in%nj)]
              }
            }else{
              overlap[[l]]=NA
            }
          }else{
            overlap[[l]]=NA
          }
        }
      }else if(all(constraints!=1)){#no connections across traps
        Gibbs="simple"
      }else{
        Gibbs="RC"
        nj=colSums(y.unk)
        trapcaps=which(nj>0)
      }
    }
    ##conditional poisson stuff  P( X = x | x >= y )
    conditionalPois=function(x,lambda,value){
      probs=rep(NA,length(x))
      for(i in 1:length(x)){
        if(value[i]>0){
          probs[i]=dpois(x[i],K*lambda[i])/(1-ppois(value[i]-1,K*lambda[i]))
        }else{
          probs[i]=dpois(x[i],K*lambda[i])
        }
      }
      return((log(probs)))
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
          if(runif(1) < exp(llycandsum-llysum)){
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
              y.true[,j]=y.true[,j]+rmultinom((n+1):M,nj[j],switchprobs)
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
                # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
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
                # switchprobs[setdiff((n+1):M,unique(ID[match]))]=0
                # switchprobs[setdiff((n+1):nUnk,unique(ID[match]))]=0  #havent' checked this yet. need to allow augmented guys
                if(length(rem)>0){
                  switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
                }
              }
              if(y.unk[l,nj]>1){
                switchprobs=switchprobs^y.unk[l,nj]
              }
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
            
            
          }else if(1==1){#complex- Fuck! ID connections across traps! 
            #calculate all switchprobs ahead of time
            switchprobs=z*lamd
            for(l in 1:nUnk){
              nj=which(y.unk[l,]>0)
              # if(length(nj)==1){
              #   switchprobs=z[(n+1):M]*K*lamd[(n+1):M,nj]
              #   match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
              #   if(length(match)!=nUnk){
              #     # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
              #     match=match[-which(ID[match]==ID[l])]#Can always match with your current cluster members
              #     legal=rep(TRUE,length(match))
              #     for(i in 1:length(match)){
              #       newcluster=which(ID==ID[match[i]])#who else is linked to the cells you are compatible with?
              #       if(length(newcluster)>0){
              #         for(j in 1:length(newcluster)){
              #           if(constraints[l,newcluster[j]]==0){#are you compatible with their cell, too?
              #             legal[i]=FALSE
              #           }
              #         }
              #       }
              #     }
              #     rem=match[!legal]
              #     match=c(l,match[legal])#take out illegal moves and add back in current spot
              #     # switchprobs[setdiff((n+1):M,unique(ID[match]))]=0
              #     # switchprobs[setdiff((n+1):nUnk,unique(ID[match]))]=0  #havent' checked this yet. need to allow augmented guys
              #     if(length(rem)>0){
              #       switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
              #     }
              #   }
              #   switchprobs=switchprobs/sum(switchprobs)
              #   if(y.unk[l,nj]>1){
              #     switchprobs=switchprobs^y.unk[l,nj]
              #     switchprobs=switchprobs/sum(switchprobs)
              #   }
              #   newID=sample((n+1):M,1,prob=switchprobs)
              # }else{
              switchprobs2=switchprobs
                match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
                if(length(match)!=nUnk){
                  # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
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
                  # switchprobs[setdiff((n+1):M,unique(ID[match]))]=0
                  # switchprobs[setdiff((n+1):nUnk,unique(ID[match]))]=0  #havent' checked this yet. need to allow augmented guys
                  if(length(rem)>0){
                    switchprobs2[rem,]=0  #havent' checked this yet. need to allow augmented guys
                  }
                }
                switchprobs2=switchprobs2/sum(switchprobs2)
                for(k in 1:length(nj)){
                  if(y.unk[l,nj[k]]>1){
                    switchprobs2[,nj[k]]=switchprobs2[,nj[k]]^y.unk[l,nj[k]]
                    # switchprobs2[,nj[k]]=switchprobs2[,nj[k]]/sum(switchprobs2[,nj[k]])
                  }
                }
                switchprobs2=switchprobs2/sum(switchprobs2)
                
                if(length(nj)==1){
                  switchprobsC=switchprobs2[,nj]/sum(switchprobs2[,nj])
                }else{
                  switchprobsC=apply(switchprobs2[,nj],1,sum)/sum(switchprobs2[,nj])
                  switchprobsC=apply(switchprobs2[,nj],1,prod)
                  switchprobsC=switchprobsC/sum(switchprobsC)
                  
                }
            
              
              # switchprobs=matrix(NA,nrow=M-n,ncol=length(nj))
              
                # for(k in 1:length(nj)){
                #   switchprobs[,k]=z[(n+1):M]*K*lamd[(n+1):M,nj[k]]
                #   match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
                #   if(length(match)!=nUnk){
                #     # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
                #     match=match[-which(ID[match]==ID[l])]#Can always match with your current cluster members
                #     legal=rep(TRUE,length(match))
                #     for(i in 1:length(match)){
                #       newcluster=which(ID==ID[match[i]])#who else is linked to the cells you are compatible with?
                #       if(length(newcluster)>0){
                #         for(j in 1:length(newcluster)){
                #           if(constraints[l,newcluster[j]]==0){#are you compatible with their cell, too?
                #             legal[i]=FALSE
                #           }
                #         }
                #       }
                #     }
                #     rem=match[!legal]
                #     match=c(l,match[legal])#take out illegal moves and add back in current spot
                #     # switchprobs[setdiff((n+1):M,unique(ID[match]))]=0
                #     # switchprobs[setdiff((n+1):nUnk,unique(ID[match]))]=0  #havent' checked this yet. need to allow augmented guys
                #     if(length(rem)>0){
                #       switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
                #     }
                #   }
                #   switchprobs[,k]=switchprobs[,k]/sum(switchprobs[,k])
                #   if(y.unk[l,nj[k]]>1){
                #     switchprobs[,k]=switchprobs[,k]^y.unk[l,nj[k]]
                #     switchprobs[,k]=switchprobs[,k]/sum(switchprobs[,k])
                #   }
                # }
                # switchprobs=apply(switchprobs,1,prod)
                # switchprobs=switchprobs/sum(switchprobs)
                
                # if(any(!is.na(overlap[[l]]))){#does this sample overlap with anyone?
                #   culprits=overlap[[l]]
                #   modify=ID[culprits]#cell probs the culprits are currently linked to. need to modify accordingly
                #   for(i in 1:length(modify)){
                #     traps=trapcaps2[[l]][[i]]
                #     #calculate switchprobs at these traps
                #     switchprobs2=matrix(NA,nrow=M-n,ncol=length(traps))
                #     for(j in 1:length(traps)){
                #       switchprobs2[,j]=z[(n+1):M]*K*lamd[(n+1):M,traps[j]]
                #     }
                #     switchprobs2=switchprobs2/sum(switchprobs2)
                #     # switchprobs2=apply(switchprobs2,1,prod)
                #     # switchprobs2=switchprobs2/sum(switchprobs2)
                #     for(j in 1:length(traps)){
                #       ncaps=y.unk[culprits[i],traps[j]]
                #       switchprobs[modify[i]]=switchprobs[modify[i]]*switchprobs2[modify[i],j]^ncaps
                #     }
                #   }
                # }
                # switchprobs=switchprobs/sum(switchprobs)
                newID=sample((n+1):M,1,prob=switchprobsC)
              
              
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
            
          }else if(1==2){#complex2
            for(l in 1:nUnk){
              nj=which(y.unk[l,]>0)
              # if(length(nj)==1){
              #   switchprobs=z[(n+1):M]*K*lamd[(n+1):M,nj]
              #   match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
              #   if(length(match)!=nUnk){
              #     # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
              #     match=match[-which(ID[match]==ID[l])]#Can always match with your current cluster members
              #     legal=rep(TRUE,length(match))
              #     for(i in 1:length(match)){
              #       newcluster=which(ID==ID[match[i]])#who else is linked to the cells you are compatible with?
              #       if(length(newcluster)>0){
              #         for(j in 1:length(newcluster)){
              #           if(constraints[l,newcluster[j]]==0){#are you compatible with their cell, too?
              #             legal[i]=FALSE
              #           }
              #         }
              #       }
              #     }
              #     rem=match[!legal]
              #     match=c(l,match[legal])#take out illegal moves and add back in current spot
              #     # switchprobs[setdiff((n+1):M,unique(ID[match]))]=0
              #     # switchprobs[setdiff((n+1):nUnk,unique(ID[match]))]=0  #havent' checked this yet. need to allow augmented guys
              #     if(length(rem)>0){
              #       switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
              #     }
              #   }
              #   switchprobs=switchprobs/sum(switchprobs)
              #   if(y.unk[l,nj]>1){
              #     switchprobs=switchprobs^y.unk[l,nj]
              #     switchprobs=switchprobs/sum(switchprobs)
              #   }
              #   newID=sample((n+1):M,1,prob=switchprobs)
              # }else{
              
              # ID probs at njs
              # if(length(nj)==1){
              #   switchprobs=z[(n+1):M]*K*lamd[(n+1):M,nj]
              # 
              # }else{
              #   switchprobs=z[(n+1):M]*K*rowSums(lamd[(n+1):M,nj])
              # }
              #   match=which(constraints[l,]==1)#you can only be moved into these cells, but they may have incompatible cells matched to them
              #   if(length(match)!=nUnk){
              #     # match=match[switchprobs[ID[match]]>0.0000001]#cut out low switchprobs guys
              #     match=match[-which(ID[match]==ID[l])]#Can always match with your current cluster members
              #     legal=rep(TRUE,length(match))
              #     for(i in 1:length(match)){
              #       newcluster=which(ID==ID[match[i]])#who else is linked to the cells you are compatible with?
              #       if(length(newcluster)>0){
              #         for(j in 1:length(newcluster)){
              #           if(constraints[l,newcluster[j]]==0){#are you compatible with their cell, too?
              #             legal[i]=FALSE
              #           }
              #         }
              #       }
              #     }
              #     rem=match[!legal]
              #     match=c(l,match[legal])#take out illegal moves and add back in current spot
              #     if(length(rem)>0){
              #       switchprobs[rem]=0  #havent' checked this yet. need to allow augmented guys
              #     }
              #   }
              #  
              #   # if(y.unk[l,nj]>1){
              #   #   switchprobs=switchprobs^y.unk[l,nj]
              #   # }
              #   switchprobs=switchprobs/sum(switchprobs)
              #   
              #   switchprobs=switchprobs^sum(y.unk[l,nj])
              #   switchprobs=switchprobs/sum(switchprobs)
              #   newID=sample((n+1):M,1,prob=switchprobs)
              
          
              
              
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
          }
            known.vector[(n+1):M]=1*(rowSums(y.true[(n+1):M,])>0)
          
          
        }else if(1==1){ #Broken MH update
          # update=sample((n+1):M,5)
          for(l in 1:nUnk){
            # for(l in update){  
            #find who you can swap IDs to
            dv<-  sqrt( (s[ID[l],1]- s[(n+1):M,1])^2 + (s[ID[l],2] - s[(n+1):M,2])^2 )
            possible<- which(dv < swap.tol)+n
            possible=possible[-which(possible==ID[l])]
            possible=possible[z[possible]==1]
            if(length(possible)>0){#Check to see if constraints match
              thiscluster=which(ID==ID[l])
              if(length(thiscluster)>1){
                thiscluster=thiscluster[-which(rowSums(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
              }else{
                thiscluster=thiscluster[-which(sum(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
              }
              if(length(thiscluster)>0){
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
              }
              if(length(possible)>0){#Can update
                jump.probability<- 1/length(possible)
                newID=ID
                # probs=(1/dv)*z
                # probs[ID[l]]=0
                # probs=probs/sum(probs)
                # newID[l]=sample(1:M,1,prob=probs)
                if(length(possible)>1){
                  newID[l] <-  sample( possible, 1)
                }else{
                  newID[l]=possible
                }
                dv<-  sqrt( (s[newID[l],1]- s[(n+1):M,1])^2 + (s[newID[l],2] - s[(n+1):M,2])^2 )
                trash<-   which(dv < swap.tol)+n
                # backprobs=(1/dv)*z
                # backprobs[newID[l]]=0
                # backprobs=backprobs/sum(backprobs)
                trash=trash[-which(trash==newID[l])]
                trash=trash[z[trash]==1]
                backcluster=which(newID==newID[l])#which guys are in the new ID cluster?
                if(length(backcluster)>1){
                  backcluster=backcluster[-which(rowSums(constraints[backcluster,])!=nUnk)]#Don't need to check if no constraints for this sc
                }else{
                  backcluster=backcluster[-which(sum(constraints[backcluster,])==nUnk)]
                }
                # Check to see if you can swap into these clusters
                if(length(backcluster)>0){
                  legal=rep(TRUE,length(trash))
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
                }
                jump.back<-  1/length(trash)
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
                
                
                ###full conditional
                nj=which(y.unk[l,]>0)
                ll.pre=ll.post=0
                y.cand2=y.true
                y.cand2[swapped,]=y.cand
                colcounts=colSums(y.unk)
                colprobs=lamd/matrix(rep(colSums(lamd),M),nrow=M,byrow=TRUE)
                # rowcounts=rowSums(y.true)
                # rowcounts.cand=rowSums(y.cand2)
                # rowprobs=rowSums(lamd)/sum(lamd)
                j=1
                ll.pre=dmultinom(y.true[,nj[j]],colcounts[nj[j]],colprobs[,nj[j]],log=TRUE)
                ll.post=dmultinom(y.cand2[,nj[j]],colcounts[nj[j]],colprobs[,nj[j]],log=TRUE)
                
                # ll.pre2=dmultinom(rowcounts,sum(rowcounts),rowprobs,log=TRUE)
                # ll.post2=dmultinom(rowcounts.cand,sum(rowcounts),rowprobs,log=TRUE)
       
           
            
                # if(rowSums(y.cand)[1]<1){#if we are emptying cluster
                  # if(rowSums(y.cand)[1]<1&rowSums(y.true[swapped,])[2]>0){#if we are emptying cluster and not creating a new one
                    
                # if(rowSums(y.true[swapped,])[2]>0){#if we are not recruiting a cluster
                if(rowSums(y.true[swapped,])[2]<(-1)){#abort conditional s update
                  #Propose new S from full conditional based on new ID
                  s.cand=matrix(NA,nrow=2,ncol=2)
                  lamd.cand=lamd
                  dtmp=matrix(NA,nrow=2,ncol=J)
                  for(i in 1:2){
                    #find mean cluster location to center window around
                    traps=which(y.cand[i,]>0)
                    if(length(traps)==0){#if empty cluster, center around current s
                      meanS=s[swapped[i],]
                    }else{#cluster not empty
                      currS=X[traps,]
                      if(!is.null(dim(currS))){#if more than one location
                        meanS=colMeans(currS)
                      }else{#only 1 location
                        meanS=currS
                      }
                    }
                    xlim2=c(meanS[1]-s.lim/2,meanS[1]+s.lim/2)
                    ylim2=c(meanS[2]-s.lim/2,meanS[2]+s.lim/2)
                    grid.locs=as.matrix(expand.grid(seq(xlim2[1],xlim2[2],res),seq(ylim2[1],ylim2[2],res)))
                    #Check for any off state space
                    rem=which(grid.locs[,1]<xlim[1]|grid.locs[,1]>xlim[2]|grid.locs[,2]<ylim[1]|grid.locs[,2]>ylim[2])
                    if(length(rem)>0){
                      grid.locs=grid.locs[-rem,]
                    }
                    dtmp.tmp=e2dist(grid.locs,X)
                    lamd.tmp=lam0*exp(-dtmp.tmp*dtmp.tmp/(2*sigma*sigma))
                    npoints=nrow(grid.locs)
                    grid.lik=rep(NA,nrow=npoints)
                    for(j in 1:npoints){
                      grid.lik[j]=exp(sum(dpois(y.cand[i,],lamd.tmp[j,],log=TRUE)))
                    }
                    grid.probs=grid.lik/sum(grid.lik)
                    choose=sample(1:npoints,1,prob=grid.probs)
                    # choose=which(grid.probs==max(grid.probs))
                    if(length(choose>1)){
                      choose=choose[1]
                    }
                    s.cand[i,]=grid.locs[choose,]
                    dtmp[i,]=dtmp.tmp[choose,]
                    lamd.cand[swapped[i],]=lamd.tmp[choose,]
                    # image(unique(grid.locs[,1]),unique(grid.locs[,2]),matrix(grid.probs,nrow=length(unique(grid.locs[,1]))))
                    # points(X)
                    # points(X[ which(y.cand[i,]>0),1],X[which(y.cand[i,]>0),2],pch=4,cex=2)
                    # points(s.cand[i,1],s.cand[i,2])
                  }
                  #update z likelihood
                  if(obstype=="poisson"){
                    pd.curr=1-exp(-lamd[swapped,])
                    pd.tmp=1-exp(-lamd.cand[swapped,])
                  }
                  pbar.tmp=(1-pd.tmp)^K
                  pbar.curr=(1-pd.curr)^K
                  prob0.tmp<- exp(rowSums(log(pbar.tmp)))
                  prob0.curr<- exp(rowSums(log(pbar.curr)))
                  fc.tmp<- prob0.tmp*psi/(prob0.tmp*psi + 1-psi)
                  fc.curr<- prob0.curr*psi/(prob0.curr*psi + 1-psi)
                  ll.z.cand=dbinom(1,1,fc.tmp)
                  ll.z.curr=dbinom(1,1,fc.curr)

                ll.y.cand=ll.y
                if(obstype=="poisson"){
                  ll.y.cand[swapped,]=dpois(y.cand,K*lamd.cand[swapped,],log=TRUE)
                }else{
                  pd.cand[swapped,]=1-exp(-lamd.cand[swapped,])
                  ll.y.cand[swapped,]=dbinom(y.cand,K,pd.cand[swapped,],log=TRUE)
                }
                

                # ll.y.cand[swapped[1],]=dpois(y.cand[1,],K*lamd.cand[swapped[1],],log=TRUE)
                # ll.y.cand[swapped[2],]=dpois(y.cand[2,],K*lamd[swapped[2],],log=TRUE)
                
                if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.z.cand))-(sum(ll.y[swapped,])+sum(ll.z.curr)))*(jump.back/jump.probability)){
                  y.true[swapped,]=y.cand
                  ll.y[swapped,]=ll.y.cand[swapped,]
                  ID[l]=newID[l]
                  known.vector[swapped]=1*(rowSums(y.cand)>0)
                  # lamd[swapped[1],]=lamd.cand[swapped[1],]
                  # s[swapped[1],]=s.cand[1,]
                  # D[swapped[1],]=dtmp[1,]
                  lamd[swapped,]=lamd.cand[swapped,]
                  s[swapped,]=s.cand
                  D[swapped,]=dtmp
                }
                
                }else{
                if(obstype=="poisson"){
                  ll.y.cand[swapped,]=dpois(y.cand,K*lamd[swapped,],log=TRUE)
                }else{
                  ll.y.cand[swapped,]=dbinom(y.cand,K,pd[swapped,],log=TRUE)
                }

                  #conditional poisson likelihood

                  # # pre
                  # ll.y.pre=conditionalPois(y.true[swapped[1],],lamd[swapped[1],],y.cand[1,])[which(y.unk[l,]>0)]
                  # #post
                  # ll.y.post=conditionalPois(y.cand[2,],lamd[swapped[2],],y.true[swapped[2],])[which(y.unk[l,]>0)]
                  # ll.y.pre[which(!is.finite(ll.y.pre))]=0
                  # ll.y.post[which(!is.finite(ll.y.post))]=0
                  # Additive poisson
                  # nj=which(y.unk[l,]>0)
                  # if(length(nj)>1){
                  #   switchprobs=z[(n+1):M]*K*rowSums(lamd[,nj])
                  #   pre=post=rowSums(y.true[,nj])
                  # }else{
                  #   switchprobs=z[(n+1):M]*K*lamd[,nj]
                  #   pre=post=y.true[,nj]
                  # }
                  # switchprobs=switchprobs/sum(switchprobs)
                  # post[ID[l]]=post[ID[l]]-sum(y.unk[l,nj])
                  # post[newID[l]]=post[newID[l]]+sum(y.unk[l,nj])
                  # pre=post=rep(0,M-n)
                  # pre[ID[l]]=sum(y.unk[l,nj])
                  # post[newID[l]]=sum(y.unk[l,nj])
                  # ll.y.pre=dmultinom(pre,sum(pre),switchprobs,log=TRUE)
                  # ll.y.post=dmultinom(post,sum(pre),switchprobs,log=TRUE)
                  # 
                  #Trap by trap MN 
                  # nj=which(y.unk[l,]>0)
                  # if(length(nj)>1){
                  #   switchprobs=pre=post=matrix(0,nrow=M-n,ncol=length(nj))
                  #   ll.y.pre=ll.y.post=rep(NA,length(nj))
                  #   for(j in 1:length(nj)){
                  #     switchprobs[,j]=z[(n+1):M]*K*lamd[,nj[j]]
                  #     switchprobs[,j]=switchprobs[,j]/sum(switchprobs[,j])
                  #     pre[ID[l],j]=y.unk[l,nj[j]]
                  #     post[newID[l],j]=y.unk[l,nj[j]]
                  #     ll.y.pre[j]=dmultinom(pre[,j],sum(pre[,j]),switchprobs[,j],log=TRUE)
                  #     ll.y.post[j]=dmultinom(post[,j],sum(pre[,j]),switchprobs[,j],log=TRUE)
                  # 
                  #   }
                  # }else{
                  #   switchprobs=z[(n+1):M]*K*lamd[,nj]
                  #   pre=post=y.true[,nj]
                  #   switchprobs=switchprobs/sum(switchprobs)
                  #   pre=post=rep(0,M-n)
                  #   pre[ID[l]]=sum(y.unk[l,nj])
                  #   post[newID[l]]=sum(y.unk[l,nj])
                  #   ll.y.pre=dmultinom(pre,sum(pre),switchprobs,log=TRUE)
                  #   ll.y.post=dmultinom(post,sum(pre),switchprobs,log=TRUE)
                  # }
                 
                  #Trap by trap MN2 (all latent n_j for focal n_j) Same LL diff as poisson
                  # nj=which(y.unk[l,]>0)
                  # if(length(nj)>1){
                  #   switchprobs=pre=post=matrix(0,nrow=M-n,ncol=length(nj))
                  #   ll.y.pre=ll.y.post=rep(NA,length(nj))
                  #   for(j in 1:length(nj)){
                  #     switchprobs[,j]=z[(n+1):M]*K*lamd[,nj[j]]
                  #     switchprobs[,j]=switchprobs[,j]/sum(switchprobs[,j])
                  #     for(k in 1:length(ID)){
                  #       pre[ID[k],j]=pre[ID[k],j]+y.unk[k,nj[j]]
                  #       post[ID[k],j]=post[ID[k],j]+y.unk[k,nj[j]]
                  #     }
                  #     post[ID[l],j]=post[ID[l],j]-y.unk[l,nj[j]]
                  #     post[newID[l],j]=post[newID[l],j]+y.unk[l,nj[j]]
                  #     ll.y.pre[j]=dmultinom(pre[,j],sum(pre[,j]),switchprobs[,j],log=TRUE)
                  #     ll.y.post[j]=dmultinom(post[,j],sum(post[,j]),switchprobs[,j],log=TRUE)
                  # 
                  #   }
                  # }else{
                  #   switchprobs=z[(n+1):M]*K*lamd[,nj]
                  #   pre=post=y.true[,nj]
                  #   switchprobs=switchprobs/sum(switchprobs)
                  #   post[ID[l]]=post[ID[l]]-sum(y.unk[l,nj])
                  #   post[newID[l]]=post[newID[l]]+sum(y.unk[l,nj])
                  #   ll.y.pre=dmultinom(pre,sum(pre),switchprobs,log=TRUE)
                  #   ll.y.post=dmultinom(post,sum(pre),switchprobs,log=TRUE)
                  # }

                  # ll.y.pre=dpois(y.unk[l,],K*lamd[swapped[1],],log=TRUE)[y.unk[l,]>0]
                  # ll.y.post=dpois(y.unk[l,],K*lamd[swapped[2],],log=TRUE)[y.unk[l,]>0]
                  
                  #IDprobs
                  # nj=which(y.unk[l,]>0)
                  # if(length(nj)>1){
                  #   switchprobs=matrix(NA,nrow=length(ID),ncol=length(nj))
                  #   for(j in 1:length(nj)){
                  #     switchprobs[,j]=lamd[ID,nj[j]]/sum(lamd[ID,nj[j]])
                  #   }
                  #   switchprobs=apply(switchprobs,1,prod)
                  #   switchprobs=switchprobs/sum(switchprobs)
                  # }else{
                  #   switchprobs=lamd[ID,nj]/sum(lamd[ID,nj])
                  # }
                  # switchprobs=lamd[ID[l],]/sum(lamd[ID[l],])
                  # switchprobs2=lamd[newID[l],]/sum(lamd[newID[l],])
                  # # 
                  # nj=which(colSums(y.cand)>0)
                  # switchprobs3=lamd/sum(lamd)
                  # switchprobsC=matrix(NA,nrow=M-n,ncol=length(nj))
                  # for(k in 1:M){
                  #   switchprobsC[k,1:length(nj)]=switchprobs3[k,nj]/sum(switchprobs3[k,nj])
                  # }
                  # switchprobsC=apply(switchprobs3[,nj],1,sum)/sum(switchprobs3[,nj])

# 
#                   ll.y.pre=dmultinom(y.unk[l,nj],sum(y.unk[l,nj]),switchprobs,log=TRUE)
#                   ll.y.post=dmultinom(y.unk[l,nj],sum(y.unk[l,nj]),switchprobs2,log=TRUE)

                  # ll.y2.pre=dmultinom(y.true[swapped[1],nj],sum(y.true[swapped[1],nj]),switchprobsC[swapped[1],],log=TRUE)
                  # ll.y2.pre2=dmultinom(y.true[swapped[2],nj],sum(y.true[swapped[2],nj]),switchprobsC[swapped[2],],log=TRUE)
                  # ll.y2.post=dmultinom(y.cand[1,nj],sum(y.cand[1,nj]),switchprobsC[swapped[1],],log=TRUE)
                  # ll.y2.post2=dmultinom(y.cand[2,nj],sum(y.cand[2,nj]),switchprobsC[swapped[2],],log=TRUE)
                # (ll.y.post+ll.y.post2)-(ll.y.pre+ll.y.pre2)
                  # switchprobs3=lamd/sum(lamd)
                  
                  # switchprobs3=lamd[swapped,]/sum(lamd[swapped,])

                  # ll.y2.pre=dmultinom(y.true[swapped[1],],sum(y.true[swapped[1],]),switchprobs3[swapped[1],],log=TRUE)
                  # ll.y2.pre2=dmultinom(y.true[swapped[2],],sum(y.true[swapped[2],]),switchprobs3[swapped[2],],log=TRUE)
                  # ll.y2.post=dmultinom(y.cand[1,],sum(y.cand[1,]),switchprobs3[1,],log=TRUE)
                  # ll.y2.post2=dmultinom(y.cand[2,],sum(y.cand[2,]),switchprobs3[2,],log=TRUE)
                  # 
                  # 
                  # #
                  # ll.y2.pre=ll.y2.pre+ll.y2.pre2
                  # ll.y2.post=ll.y2.post+ll.y2.post2
                  
                  # if(runif(1)<exp((sum(ll.y.post)+sum(ll.y2.post))-(sum(ll.y.pre)+sum(ll.y2.pre)))*(jump.back/jump.probability)){
                    
                if(runif(1)<exp(sum(ll.y2.post)-sum(ll.y2.pre))*(jump.back/jump.probability)){
                  # if(runif(1)<exp((sum(ll.y.post)+sum(ll.y.cand[swapped,]))-(sum(ll.y.pre)+sum(ll.y[swapped,])))*(jump.back/jump.probability)){
                    
                    # if(runif(1)<exp((sum(ll.y.cand[swapped,]))-(sum(ll.y[swapped,])))*(jump.back/jump.probability)){
                  # if(runif(1)<exp((sum(ll.y.cand[swapped,]))-(sum(ll.y[swapped,])))*(backprobs[swapped[1]]/probs[swapped[2]])){
                  # if(runif(1)<exp((sum(ll.y.cand))-(sum(ll.y)))*(jump.back/jump.probability)){
                # if(runif(1)<exp((sum(ll.y.cand[swapped,])+sum(ll.z.cand))-(sum(ll.y[swapped,])+sum(ll.z.curr)))*(jump.back/jump.probability)){
                  y.true[swapped,]=y.cand
                  ll.y[swapped,]=ll.y.cand[swapped,]
                  ID[l]=newID[l]
                  known.vector[swapped]=1*(rowSums(y.cand)>0)
                  # ll.y=ll.y.cand
                  # lamd=lamd.cand
                  # lam0=lam0.cand
                  # sigma=sigma.cand
                  # Update lam0
                  # for(up in 1:25){
                  #   llysum=sum(ll.y)
                  #   lam0.cand<- rnorm(1,lam0,proppars$lam0/2)
                  #   if(lam0.cand > 0){
                  #     lamd.cand<- lam0.cand*exp(-D*D/(2*sigma*sigma))
                  #     ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
                  #     llycandsum=sum(ll.y.cand)
                  #     if(runif(1) < exp(llycandsum-llysum)){
                  #       lam0<- lam0.cand
                  #       lamd=lamd.cand
                  #       ll.y=ll.y.cand
                  #     }
                  #   }
                  #   for (i in 1:2) {
                  #     Scand <- c(rnorm(1, s[swapped[i], 1], proppars$sx), rnorm(1, s[swapped[i], 2], proppars$sy))
                  #     if(useverts==FALSE){
                  #       inbox <- Scand[1] < xlim[2] & Scand[1] > xlim[1] & Scand[2] < ylim[2] & Scand[2] > ylim[1]
                  #     }else{
                  #       inbox=inout(Scand,vertices)
                  #     }
                  #     if (inbox) {
                  #       dtmp <- sqrt((Scand[1] - X[, 1])^2 + (Scand[2] - X[, 2])^2)
                  #       lamd.cand[swapped[i],]<- lam0*exp(-dtmp*dtmp/(2*sigma*sigma))
                  #       if(obstype=="bernoulli"){
                  #         pd.cand[swapped[i],]=1-exp(-lamd.cand[swapped[i],])
                  #         ll.y.cand[swapped[i],]= dbinom(y.true[swapped[i],],K,pd.cand[swapped[i],]*z[swapped[i]],log=TRUE)
                  #         if (runif(1) < exp(sum(ll.y.cand[swapped[i],]) - sum(ll.y[swapped[i],]))) {
                  #           s[swapped[i],]=Scand
                  #           D[swapped[i],]=dtmp
                  #           lamd[swapped[i],]=lamd.cand[swapped[i],]
                  #           pd[swapped[i],]=pd.cand[swapped[i],]
                  #           ll.y[swapped[i],]=ll.y.cand[swapped[i],]
                  #         }
                  #       }else{#poisson
                  #         ll.y.cand[swapped[i],]= dpois(y.true[swapped[i],],K*lamd.cand[swapped[i],]*z[swapped[i]],log=TRUE)
                  #         if (runif(1) < exp(sum(ll.y.cand[swapped[i],]) - sum(ll.y[swapped[i],]))) {
                  #           s[swapped[i],]=Scand
                  #           D[swapped[i],]=dtmp
                  #           lamd[swapped[i],]=lamd.cand[swapped[i],]
                  #           ll.y[swapped[i],]=ll.y.cand[swapped[i],]
                  #         }
                  #       }
                  #     }
                  #   }
                  # sigma.cand<- rnorm(1,sigma,proppars$sigma/2)
                  # if(sigma.cand > 0){
                  #   lamd.cand<- lam0*exp(-D*D/(2*sigma.cand*sigma.cand))
                  #   ll.y.cand= dpois(y.true,K*lamd.cand*z,log=TRUE)
                  #   llycandsum=sum(ll.y.cand)
                  #   if(runif(1) < exp(llycandsum-llysum)){
                  #     sigma<- sigma.cand
                  #     lamd=lamd.cand
                  #     ll.y=ll.y.cand
                  #   }
                  # }
                  # }
                }
                }
              }
            }
          }
        }else{ #weird shit
        for(l in 1:nUnk){
          #find who you can swap IDs to
          possible=(n+1):M
          possible=possible[z[possible]==1]
          if(length(possible)>0){#Check to see if constraints match
            thiscluster=which(ID==ID[l])
            if(length(thiscluster)>1){
              thiscluster=thiscluster[-which(rowSums(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
            }else{
              thiscluster=thiscluster[-which(sum(constraints[thiscluster,])==nUnk)]#Don't need to check if no constraints for this sc
            }
            if(length(thiscluster)>0){
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
            }
            if(length(possible)>0){#Can update
              newID=ID
              ll.y.tmp=rep(-Inf,M-n)
              y.cand=matrix(0,nrow=M-n,ncol=J)
              
              for(l2 in possible){
                newID[l]=l2
                clusterguys=which(newID==l2)
                if(length(clusterguys)>0){
                  if(length(clusterguys)>1){
                    y.cand[l2,]=colSums(y.unk[clusterguys,])
                  }else{
                    y.cand[l2,]=y.unk[clusterguys,]
                  }
                }
                ll.y.tmp[l2]=sum(dpois(y.cand[l2,],K*lamd[l2,],log=TRUE))
              }
              swap.probs=exp(ll.y.tmp)/sum(exp(ll.y.tmp))
              newID[l]=sample((n+1):M,1,prob=swap.probs)
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
            
            y.true[swapped,]=y.cand
            ll.y[swapped,]=ll.y.cand[swapped,]
            ID[l]=newID[l]
            known.vector[swapped]=1*(rowSums(y.cand)>0)
            }
          }
        }
        }
      }
      
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
        out[idx,]<- c(lam0,sigma ,sum(z),length(unique(ID)))
        # print(out[idx,])
        # out[idx,]<- c(lam0,sigma ,sum(z), sum(rowSums(y.true[(n+1):M,])>0))
        idx=idx+1
      }
    }  # end of MCMC algorithm

    if(keepACs==TRUE){
      list(out=out, sxout=sxout, syout=syout, zout=zout,IDout=IDout)
    }else{ 
      list(out=out)
    }
  }





