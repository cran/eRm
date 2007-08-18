print.threshold <- function(x,...)
{
    mt_vek <- apply(x$X, 2, max, na.rm = TRUE)
    categ <- sequence(mt_vek)
    it <- rep(1:ncol(x$X), mt_vek)
    categ.l <-split(x$threshpar,it)

    thrpar <- t(sapply(categ.l,"[",1:max(categ)))
    itloc<-apply(thrpar,1,mean,na.rm=TRUE)
    thrcoef<-as.data.frame(cbind(colnames(x$X),round(itloc,5),round(thrpar,5)))
    colnames(thrcoef)<-c("Item","Location",paste("Threshold",1:max(categ),sep=""))
    cat("\n")
    print(thrcoef)
    cat("\n")
}