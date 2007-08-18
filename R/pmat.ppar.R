`pmat.ppar` <-
function(object)
# computes a list of expected probabilities for objects of class "ppar" for each NA-subgroup
# without category!
{

X <- object$X
mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
mt_ind <- rep(1:length(mt_vek),mt_vek)

rp <- rowSums(X,na.rm=TRUE)
maxrp <- sum(mt_vek)
TFrow <- ((rp==maxrp) | (rp==0))
#if (any(TFrow)) {
  #cat("Warning message: For the following persons no expected probabilites are computed due to 0/full raw score: \n")
#  cat(rownames(object$X)[TFrow],sep=", ")
#  cat("\n")
  #X <- X[!TFrow,]
#}



pmat.l <- lapply(object$thetapar, function(theta1) {   
             theta <- theta1
             p.list <- tapply(object$betapar,mt_ind,function(beta.i) {     #matrices of expected prob as list (over items)
                     beta.i <- c(0,beta.i)
                     ind.h <- 0:(length(beta.i)-1)
                     theta.h <- ind.h %*% t(theta)
                     tb <- exp(theta.h+beta.i)
                     denom <- colSums(tb)
                     pi.mat <- apply(tb,1,function(y) {y/denom})
                     return(pi.mat)
                   })
    p.list0 <- lapply(p.list,function(pl) {pl[,-1]})               #delete 0th category
    pmat <- matrix(unlist(p.list0),nrow=length(theta1))      #save as matrix
    return(pmat)
}) 

#----------labels----------
#names(pmat.l) <- paste("NAgroup",1:length(pmat.l),sep="")
cnames <- substr(names(object$betapar),6,40)
for (i in 1:length(pmat.l)) {
      dimnames(pmat.l[[i]]) <- list(names(object$thetapar[[i]]),cnames)}
#-----------end labels-------   
NApos <- tapply(1:length(object$gmemb),object$gmemb,function(ind) {   #positions for NA replacement
                       xvec <- X[ind,][1,]
                       which(is.na(xvec))
                       })
pmat <- NULL
for (i in 1:length(pmat.l)) {
       pmat.l[[i]][,NApos[[i]]] <- NA            #insert NA's
       pmat <- rbind(pmat,pmat.l[[i]])
       }
return(pmat)
}

