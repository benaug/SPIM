LRmatchOpen <-  function(M, left, nleft, right,nright, X, Nfixed){
  # This function takes a left and right data set and ad hocly associates unknown right guys
  # with unknown left guys to minimize the total separation distance.
  #Unknown left and right only guys can't be captured both guys because then they would be known.
  #Function modified from previous version by Andy Royle

  # function needs to spit out initial ID and activity centers
  #put in flag for cases where left or right are all known
  idkpl<- (Nfixed+1):dim(left)[1]
  idkpr<- (Nfixed+1):dim(right)[1]
  #Extract unknown individuals
  ld <- left[idkpl,,,]
  rd <-right[idkpr,,,]
  #Let's take the year dimension and extend the traps dimension
  ld2=ld[,,,1]
  rd2=rd[,,,1]
  X2=X[[1]][,1:2]
  t=dim(left)[4]
  for(l in 2:t){
    ld2=abind(ld2,ld[,,,l],along=2)
    rd2=abind(rd2,rd[,,,l],along=2)
    X2=rbind(X2,X[[l]][,1:2])
  }
  X2=as.matrix(X2)
  #The rest of the closed pop code should now work

  #Sum counts for each unknown individual/trap
  ld2<- apply(ld2, c(1,2), sum) # total counts per trap
  rd2<- apply(rd2,c(1,2), sum)

  #matrices to store initial activity centers
  sbar.left<- matrix(NA,nrow=nleft,ncol=2)
  sbar.right<- matrix(NA,nrow=nright,ncol=2)

  #record average location of capture
  for(i in 1:nleft){
    if(sum(ld2[i,])>0){  # should always be satisfied... unless you didn't remove uncaptured individuals from simulated data set.
      pos.traps<- (1:nrow(X2))[ld2[i,]>0] #Where was this guy left captured?
      pos.traps<- rep(pos.traps, ld2[i,][ld2[i,]>0]) #rep each trap by number of captures
      if(length(pos.traps)>1){
        sbar.left[i,]<- apply(matrix(X2[pos.traps,],ncol=2,byrow=FALSE),2,mean) #record average location of capture
      }else{
        sbar.left[i,]=X2[pos.traps,]
      }
    }else{
      sbar.left[i,]<- X2[sample(1:nrow(X2),1),] #if not captured, pick random location.
    }
  }
  for(i in 1:nright){
    if(sum(rd2[i,])>0){
      pos.traps<- (1:nrow(X2))[rd2[i,]>0]
      pos.traps<- rep(pos.traps, rd2[i,][rd2[i,]>0])
      if(length(pos.traps)>1){
        sbar.right[i,]<- apply(matrix(X2[pos.traps,],ncol=2,byrow=FALSE),2,mean)
      }else{
        sbar.right[i,]=X2[pos.traps,]
      }
    }else{
      sbar.right[i,]<- X2[sample(1:nrow(X2),1),]
    }
  }

  D<- e2dist(sbar.right,sbar.left)
  # optimization problem here is to put a 1 in each row such that sum of all distance is small
  ID_R2L<- sample(1:nleft, nright)
  Q<- sum( D[cbind(1:nright,ID_R2L)]  )

  if(nleft > nright){ #Should always happen
    for(loop in 1:20){
      for(i in 1:nrow(D)){
        # if there are unused left guys then try to make a swap there first
        notused<- (1:nleft)[is.na(match(1:nleft,ID_R2L))]
        curr.spot<- ID_R2L[i]
        Qtmp<- rep(NA,length(notused))
        for(k in 1:length(notused)){
          ID_R2L[i]<- notused[k]
          Qtmp[k]<- sum(  D[cbind(1:nright,ID_R2L)] )
        }
        if(min(Qtmp) < Q ){
          # Make the swap
          swap.in<- Qtmp==min(Qtmp)
          ID_R2L[i]<- notused[Qtmp==min(Qtmp)][1]  # Just use the first one
          Q<- min(Qtmp)
        }
        else{
          ID_R2L[i]<- curr.spot
        }
      }
    }
  }
  ## for the last loop no other point could change with the available points.
  for(loop in 1:20){
    for(i in 1:nrow(D)){
      curr.spot<- ID_R2L[i]
      Qtmp<- rep(NA,length(ID_R2L))  # loop over EACH other
      for(k in 1:length(ID_R2L)){
        ID_R2L[i]<-  ID_R2L[k]
        ID_R2L[k]<- curr.spot
        Qtmp[k]<- sum(  D[cbind(1:nright,ID_R2L)] )
        ID_R2L[k]<- ID_R2L[i] # set it back to where it was after computing criterion
        ID_R2L[i]<- curr.spot
      }
      if(min(Qtmp) < Q ){
        # Make the swap
        which<- (1:length(ID_R2L))[Qtmp==min(Qtmp)][1]
        swap.in<- ID_R2L[which]
        ID_R2L[i]<- swap.in
        ID_R2L[which]<- curr.spot
        Q<- min(Qtmp)
      }else{  # Make sure to put it back if no swap was made
        ID_R2L[i]<- curr.spot
      }
      #cat("new Q: ", Q, fill=TRUE)
    }
  }
  #Randomly match lefts to boths not captured
  ID_L=sample((Nfixed+1):M,nleft)
  #Translate new left IDs to right IDs
  ID_R=ID_L[ID_R2L]
  #Add back to known individuals
  if(Nfixed>0){
    ID_L<- c((1:Nfixed) , ID_L )
    ID_R<- c((1:Nfixed) , ID_R )
  }
  return(list(ID_L=ID_L,ID_R=ID_R))
}
