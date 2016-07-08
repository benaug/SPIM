#' Build a data set for the spatial partial ID camera trap model
#' @param input a data frame with columns ID, trap, occ, and type.  Each capture event has its own row with individual ID, capture
#' occasion, trap number, and capture type (B, L, or R)
#' @param K the integer number of capture occasions
#' @param X the J x 3 matrix of trap locations and number of cameras at each trap. columns are X, Y, #cams (1 or 2)
#' @param IDknown the vector listing the complete indentity individuals, 1:nC.  Leave blank or set to NA if no identities are complete
#' @param buff the distance to buffer the trapping array in the X and Y dimensions to produce the state space
#' @param vertices a matrix of n_verts X 2 X and Y locations for the vertices of a polygon state space
#' @author Ben Augustine
#' @description This function formats the input object into a data set in the necessary format to run the spatial partial ID model.  You should
#' number the nC complete identity individuals as 1:nC with individuals that had a both capture numbered 1:nB and any other complete identity
#' individuals numbered (nB+1):nC.
#' @return a list with elements pIDL and pIDR, which are also lists.  Element X of pIDL corresponds to left individual X and contains the right
#' individuals it was matched with on at least one MCMC iteration and the posterior for each match.  Similarly for pIDR.
#' @examples
#' \dontrun{
#' #trivial example
#' X=cbind(1:6,1:6,rep(2,6))
#' ID=c(1,1,2,4,5,3)
#' occ=c(1,2,1,4,2,1)
#' trap=c(1,2,3,3,2,1)
#' type=c("B","L","B","R","R","L")
#' input=data.frame(ID=ID,trap=trap,occ=occ,type=type)
#' data=build.data(input,X=X,K=5,IDknown=1:2,buff=2,model="2side")
#'}

build.data=function(input,K,X,IDknown=NA,buff=NA,vertices=NA,model="2side"){
  #Check to see if IDs ordered correctly
  IDs=sort(unique(input$ID))
  if(all(IDs!=1:length(IDs))){
    stop("Individuals not ordered consecutively starting at 1")
  }
  n=max(input[,1])
  J=nrow(X)
  if(model=="2side"){
    B=input[input[,4]=="B",]
    L=input[input[,4]=="L",]
    R=input[input[,4]=="R",]
    nB=nrow(B)
    nL=nrow(L)
    nR=nrow(R)
    Bids=sort(unique(B$ID))
    nb=length(Bids)
    if(all(is.na(IDknown))){
      IDknown=integer(0)
    }
    nC=length(IDknown)#number of complete identity guys
    #Make sure B guys are in IDknown and are numbered 1:nb and known L or R guys are numbered nb+1, nb+2, ...
    if(nB>0){#if we have a both data set
      if(any(Bids!=1:nb)){
        stop("ID numbers for both captures must match 1:nb entries in IDknown and IDknown must be 1:nb")
      }
      if(nC>nb){#if any of the L or R only guys are known
        nremain=nC-nb
        if(any(IDknown[(nb+1):(nb+nremain)]!=(nb+1):(nb+nremain))){
          stop("ID numbers for any complete identity individuals must be numbered nb+1, nb+2, ... ")
        }
      }
    }
    if(nB==0&nC>0){#if we have any complete identity individuals but no B captures
      if(any(IDknown!=1:nC)){
        stop("ID numbers for complete identity individuals must match the 1:nC entries in IDknown and IDknown must be 1:nC")
      }
    }
    #Make data sets
    both=array(0,dim=c(nB,3,K,J))
    left=array(0,dim=c(nL,3,K,J))
    right=array(0,dim=c(nR,3,K,J))
    if(nB>0){
      for(i in 1:nB){
        both[B[i,1],1,B[i,3],B[i,2]]=1
      }
    }
    if(nL>0){
      for(i in 1:nL){
        left[L[i,1],2,L[i,3],L[i,2]]=1
      }
    }
    if(nR>0){
      for(i in 1:nR){
        right[R[i,1],3,R[i,3],R[i,2]]=1
      }
    }
    both=both[rowSums(both)>0,,,]
    nb=nrow(both)
    #Remove 0s from left and right
    idx=which(rowSums(left)==0)
    idx=idx[idx>nb]
    if(length(idx)>0){
      left=left[-idx,,,]
    }
    idx=which(rowSums(right)==0)
    idx=idx[idx>nb]
    if(length(idx)>0){
      right=right[-idx,,,]
    }
    #IDs
    ID_L=rep(NA,nrow(left))
    ID_R=rep(NA,nrow(right))
    if(nC>0){
      ID_L[1:nC]=ID_R[1:nC]=1:nC
    }
    #buffer or vertices
    if(all(!is.na(vertices))){
      data=list(both=both,left=left,right=right,X=X,IDknown=IDknown,n=n,ID_L=ID_L,ID_R=ID_R,K=K,vertices=vertices)
    }else if (!is.na(buff)){
      data=list(both=both,left=left,right=right,X=X,IDknown=IDknown,n=n,ID_L=ID_L,ID_R=ID_R,K=K,buff=buff)
    }else{
      stop("User must input a buffer or polygon vertices")
    }
  }else if(model=="SCR"){
    y=array(0,dim=c(n,K,J))
    for(i in 1:nrow(input)){
      y[input[i,1],input[i,3],input[i,2]]=1
    }
    if(!is.na(vertices)){
      data=list(y=y,X=X,n=n,K=K,vertices=vertices)
    }else if (!is.na(buff)){
      data=list(y=y,X=X,n=n,K=K,buff=buff)
    }else{
      stop("User must input a buffer or polygon vertices")
    }
  }else{
    stop("model must be 2side or SCR")
  }
  return(data)
}
