thresholds.eRm <- function(object)                # uses matrix approach
{
#Computation of threshold parameters for polytomous models
#object of class "eRm" (but not "dRm")

  if ((object$model == "LLTM") || (object$model == "RM")) stop("Threshold parameters are computed only for polytomous models!")
  if ((object$model == "LRSM") || (object$model == "LPCM")) {
    mpoints <- object$mpoints
    ngroups <- object$ngroups
    vecrep <- mpoints * ngroups                      
  } else {
    mpoints <- 1
    ngroups <- 1
    vecrep <- 1
  }
  
  betapar <- object$betapar
  indmt <- apply(object$X,2,max,na.rm=TRUE)         #number of categories per item
  mt_vek1 <- sequence(indmt[1:(length(indmt)/mpoints)]) #1 block of beta-items
  mt_vek <- rep(mt_vek1, vecrep) 
  sq<-ifelse(mt_vek > 1,-1,0)
  d1<-diag(sq[-1])
  k<-length(betapar)
  d2<-diag(k)
  d2[-k,-1]<-d2[-k,-1]+d1
  threshpar <-as.vector(crossprod(betapar,d2)*-1)                  #vector with threshold parameters
  
  names(threshpar) <- paste("thresh",names(betapar))
  
  vc.beta <- (object$W%*%solve(object$hessian)%*%t(object$W)) #VC matrix beta's
  se.thresh <- sqrt(diag(d2%*%(vc.beta)%*%t(d2)))             #standard errors of thresholds
  names(se.thresh) <- names(threshpar)

  blocks <- rep(1:vecrep, each = length(mt_vek1))
  thblock <- split(threshpar,blocks)                          #block of threshholds (as in design matrix)
  indmt1 <- indmt[1:(length(indmt)/mpoints)]
  indvec <- rep(1:length(indmt1),indmt1)
  
  threshtab.l <- lapply(thblock, function(x) {                     #list of table-blocks
                     Location <- tapply(x,indvec,mean)             #location parameters
                     thresh.l <- split(x, indvec)
                     threshmat <- t(sapply(thresh.l,"[",1:max(mt_vek)))
                     colnames(threshmat) <- paste("Threshold", 1:dim(threshmat)[2])
                     parmat <- cbind(Location,threshmat)
                     }) 
  
  #determine item names for block-table
  cnames <- colnames(object$X)
  ind.it <- rep(1:mpoints,each = length(cnames)/mpoints)           #item label index
  itnames1 <- as.vector(unlist(tapply(cnames, ind.it, function(x) rep(x, ngroups)))) 
  rep.ind <- sapply(threshtab.l, function(x) dim(x)[1])
  sp.ind <- rep(1:length(rep.ind), rep.ind)

  names.l <- split(itnames1, sp.ind)                   #names as list
  for (i in 1:length(threshtab.l)) rownames(threshtab.l[[i]]) <- names.l[[i]]              #name the items

  result <- list(threshpar = threshpar,se.thresh = se.thresh, threshtable = threshtab.l)
  class(result) <- "threshold"
  result

}