#' Plot trapping array and realized capture results of simulated data set
#' @param data a list produced by sim2side or in the same format
#' @param offset the amount to offset each camera for double traps
#' @param cexX the cex for representing the cameras
#' @param cexB the cex for complete identity activity centers
#' @param cexS the cex for partial identity activity centers
#' @param cex0 the cex for uncaptured activity centers
#' @param cexT the cex for the text indicating B, L, R, or LR captures for each activity center
#' @param main the title for the plot
#' @author Ben Augustine
#' @description This function plots the realized capture results from the simulation function. It plots the juxtaposition of
#' traps and activity centers, differentiating between complete identity individuals (B), partial identity individuals captured
#' on the left side only (L), right side only (R), left and right side but not simultanously (LR), and uncaptured individuals.
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
#'sim2side.plot(data)
#'}

sim2side.plot=function(data,offset=0.15,cexX=1.5,cexB=3,cexS=3,cex0=1,cexT=1,main="",plotACs=TRUE,plottimes=1){
  if(any(plottimes!=1)&!("tf"%in%names(data))){
    stop("data must have element 'tf' if plottimes!=1")
  }
  if(all(plottimes==1)){#If static grid
    if("grid"%in%names(data)){ #if discrete state space
      grid=data$grid
      X=as.matrix(data$X)
      s=data$s[,2:3]
      grid2=xtabs(grid$cov~grid$x+grid$y)
      image(x=as.numeric(rownames(grid2)),y=as.numeric(colnames(grid2)),z=grid2,xlab="X",ylab="Y",main=main)
      points(X[X[,3]==1,1:2],pch=4,cex=cexX,lwd=2)
      points(X[X[,3]==2,1]-offset,X[X[,3]==2,2],pch=4,cex=cexX,lwd=2)
      points(X[X[,3]==2,1]+offset,X[X[,3]==2,2],pch=4,cex=cexX,lwd=2)
    }else{ # continuous state space
      if("vertices"%in%names(data)){stop("This function only works for rectangular continuous state spaces or discrete state spaces")}
      X=as.matrix(data$X)
      s=data$s
      xlim=range(X[,1])+c(-buff,buff)
      ylim=range(X[,2])+c(-buff,buff)
      plot(X[X[,3]==1,1:2],xlim=xlim,ylim=ylim,pch=4,cex=cexX,lwd=2,xlab="X",ylab="Y",yaxt="n",main=main)
      points(X[X[,3]==2,1]-offset,X[X[,3]==2,2],pch=4,cex=cexX,lwd=2)
      points(X[X[,3]==2,1]+offset,X[X[,3]==2,2],pch=4,cex=cexX,lwd=2)
    }

    if(plotACs==TRUE){
      points(s[data$IDknown,],col=rgb(69,139,0,200,maxColorValue = 255),pch=20,cex=cexB)
      singles=setdiff(unique(c(data$ID_L,data$ID_R)),1:nrow(data$both))
      if(length(data$IDknown)>0){
        B=max(data$IDknown)
      }else{
        B=0
      }
      nocaps=setdiff((max(B)+1):N,singles)
      points(s[nocaps,],col="black",pch=20,cex=cex0)
      points(s[singles,],col=rgb(255,195,15,200,maxColorValue = 255),pch=20,cex=cexS)
      lefts=setdiff(data$ID_L,data$ID_R)
      rights=setdiff(data$ID_R,data$ID_L)
      LRs=setdiff(base::intersect(data$ID_R,data$ID_L),1:nrow(data$both))
      if(length(lefts)>0){
        text(x=s[lefts,1],y=s[lefts,2],rep("L",length(lefts)),cex=cexT,font=2)
      }
      if(length(rights)>0){
        text(x=s[rights,1],y=s[rights,2],rep("R",length(rights)),cex=cexT,font=2)
      }
      if(length(LRs)){
        text(x=s[LRs,1],y=s[LRs,2],rep("LR",length(LRs)),cex=cexT,font=2)
      }
      if(B>0){
        text(x=s[1:B,1],y=s[1:B,2],rep("B",length(data$IDknown)),cex=cexT,font=2)
      }
    }
  }else{#If dynamic grid
    for(i in 1:length(plottimes)){
      if("grid"%in%names(data)){ #if discrete state space
        grid=data$grid
        X=as.matrix(data$X)
        s=data$s[,2:3]
        grid2=xtabs(grid$cov~grid$x+grid$y)
        image(x=as.numeric(rownames(grid2)),y=as.numeric(colnames(grid2)),z=grid2,xlab="X",ylab="Y",main=paste("k =",plottimes[i]))
        points(X[tf[,plottimes[i]]==1,1:2],pch=4,cex=cexX,lwd=2)
        points(X[tf[,plottimes[i]]==2,1]-offset,X[tf[,plottimes[i]]==2,2],pch=4,cex=cexX,lwd=2)
        points(X[tf[,plottimes[i]]==2,1]+offset,X[tf[,plottimes[i]]==2,2],pch=4,cex=cexX,lwd=2)
      }else{ # continuous state space
        if("vertices"%in%names(data)){stop("This function only works for rectangular continuous state spaces or discrete state spaces")}
        X=as.matrix(data$X)
        s=data$s
        xlim=range(X[,1])+c(-buff,buff)
        ylim=range(X[,2])+c(-buff,buff)
        plot(X[tf[,plottimes[i]]==1,1:2],xlim=xlim,ylim=ylim,pch=4,cex=cexX,lwd=2,xlab="X",ylab="Y",yaxt="n",main=paste("k =",plottimes[i]))
        points(X[tf[,plottimes[i]]==2,1]-offset,X[tf[,plottimes[i]]==2,2],pch=4,cex=cexX,lwd=2)
        points(X[tf[,plottimes[i]]==2,1]+offset,X[tf[,plottimes[i]]==2,2],pch=4,cex=cexX,lwd=2)
      }

      if(plotACs==TRUE){
        points(s[data$IDknown,],col=rgb(69,139,0,200,maxColorValue = 255),pch=20,cex=cexB)
        singles=setdiff(unique(c(data$ID_L,data$ID_R)),1:nrow(data$both))
        if(length(data$IDknown)>0){
          B=max(data$IDknown)
        }else{
          B=0
        }
        nocaps=setdiff((max(B)+1):N,singles)
        points(s[nocaps,],col="black",pch=20,cex=cex0)
        points(s[singles,],col=rgb(255,195,15,200,maxColorValue = 255),pch=20,cex=cexS)
        lefts=setdiff(data$ID_L,data$ID_R)
        rights=setdiff(data$ID_R,data$ID_L)
        LRs=setdiff(base::intersect(data$ID_R,data$ID_L),1:nrow(data$both))
        if(length(lefts)>0){
          text(x=s[lefts,1],y=s[lefts,2],rep("L",length(lefts)),cex=cexT,font=2)
        }
        if(length(rights)>0){
          text(x=s[rights,1],y=s[rights,2],rep("R",length(rights)),cex=cexT,font=2)
        }
        if(length(LRs)){
          text(x=s[LRs,1],y=s[LRs,2],rep("LR",length(LRs)),cex=cexT,font=2)
        }
        if(B>0){
          text(x=s[1:B,1],y=s[1:B,2],rep("B",length(data$IDknown)),cex=cexT,font=2)
        }
      }
      Sys.sleep(2)
    }

  }
}
