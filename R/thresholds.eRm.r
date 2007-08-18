thresholds.eRm <- function(object)                # uses matrix approach
{
#Computation of threshold parameters for polytomous models
#object of class "eRm" (but not "dRm")

  if ((res$model == "LLTM") || (res$model == "RM")) stop("Threshold parameters are computed only for polytomous models!")

  bp<-object$betapar
  indmt <- apply(object$X,2,max,na.rm=TRUE)
  mt_vek <- sequence(indmt)
  sq<-ifelse(mt_vek > 1,-1,0)
  d1<-diag(sq[-1])
  k<-length(bp)
  d2<-diag(k)
  d2[-k,-1]<-d2[-k,-1]+d1
  threshpar <-as.vector(crossprod(bp,d2)*-1)                  #vector with threshold parameters
  names(threshpar) <- paste(rep(colnames(object$X),indmt),paste("thresh",mt_vek,sep=""))

  vc.beta <- (object$W%*%solve(object$hessian)%*%t(object$W)) #VC matrix beta's
  se.thresh <- sqrt(diag(d2%*%(vc.beta)%*%t(d2)))             #standard errors of thresholds
  names(se.thresh) <- names(threshpar)

  result <- list(threshpar = threshpar,se.thresh = se.thresh, X = object$X)
  class(result) <- "threshold"
  result

}