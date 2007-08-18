`plist.internal` <-
function(object,theta)
# computes a list of expected probabilities for objects of class Rm
# with 0th category included!
{

X <- object$X
mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
mt_ind <- rep(1:length(mt_vek),mt_vek)


#--------compute list matrix of probabilites for fixed theta) 
p.list <- tapply(object$betapar,mt_ind,function(beta.i) {
                 beta.i <- c(0,beta.i)
                 ind.h <- 0:(length(beta.i)-1)
                 theta.h <- ind.h %*% t(theta)          #multiply category with 
                 #tb <- exp(theta.h-beta.i)
                 tb <- exp(theta.h+beta.i)
                 denom <- colSums(tb)
                 pi.mat <- apply(tb,1,function(y) {y/denom})
                 return(pi.mat)
               })
return(p.list)
}

