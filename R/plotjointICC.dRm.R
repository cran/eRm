`plotjointICC.dRm` <-
function(object, item.subset = "all", legend=TRUE, xlim=c(-4,4),ylim=c(0,1),
         xlab="Latent Dimension",ylab="Probability to Solve",lty=1,...)

#produces one common ICC plot for Rasch models only
#object of class "dRm"
#item.subset...specify items that have to be plotted; if NA, all items are used
#legend...if legend should be plotted

{
  theta <- seq(xlim[1],xlim[2],by=0.1)

  if (any(item.subset=="all")) {
    it.legend <- 1:dim(object$X)[2]
  } else {
    if (is.character(item.subset)) { 
      it.legend <- item.subset
      betatemp <- t(as.matrix(object$betapar))
      colnames(betatemp) <- colnames(object$X)
      object$betapar <- betatemp[,item.subset]
    } else {
      it.legend <- colnames(object$X)[item.subset]
      object$betapar <- object$betapar[item.subset]
    }
    object$X <- object$X[,item.subset]                            #pick out items defined in itemvec 
  }


  th.ord <- order(theta)

  p.list <- plist.internal(object,theta)
  p.list <- lapply(p.list,function(x) {x[,-1]})               #Delete 0-probabilites
  p.mat <- matrix(unlist(p.list),ncol=length(p.list))
  text.ylab <- p.mat[(1:length(theta))[theta==median(theta)],]
  
  get(getOption("device"))()

  matplot(sort(theta),p.mat[th.ord,],type="l",lty=lty,col=1:(dim(p.mat)[2]),
          main=("ICC plot"),xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  if (!legend) {
    text(x=median(theta),y=text.ylab,labels=paste("I",1:(dim(p.mat)[2]),sep=""),col=1:(dim(p.mat)[2]))
  } else {
    legend(xlim[1],ylim[2],paste("Item",it.legend),col=1:(dim(p.mat)[2]),lty=lty,...)
  }
}

