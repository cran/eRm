`Waldtest.Rm` <-
function(object, splitcr="median")
{
# performs item-based Wald test (Fischer & Molenaar, p.90)
# object... object of class RM
# splitcr... splitting criterion for LR-groups. "median" to a median raw score split,
#            "mean" corobjectponds to the mean raw score split.
#            optionally also a vector of length n for group split can be submitted.


X.original<-object$X
numsplit<-is.numeric(splitcr)
if (any(is.na(object$X))) {
  if (!numsplit && splitcr=="mean") {                                   #mean split
    X<-object$X
    # calculates index for NA groups
    # from person.parameter.eRm
      dichX <- ifelse(is.na(X),1,0)
      strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
      gmemb <- as.vector(data.matrix(data.frame(strdata)))
    gindx<-unique(gmemb)
    rsum.all<-rowSums(X,na.rm=T)
    grmeans<-tapply(rsum.all,gmemb,mean)      #sorted
    ngr<-table(gmemb)                         #sorted
    m.all<-rep(grmeans,ngr)                   #sorted,expanded
    rsum.all<-rsum.all[order(gmemb)]
    spl<-ifelse(rsum.all<m.all,1,2)
    splitcr<-spl
    object$X<-X[order(gmemb),]
  }
  if (!numsplit && splitcr=="median") {                                   #median split
    cat("Warning message: Persons with median raw scores are assigned to the lower raw score group!\n")
    X<-object$X
    # calculates index for NA groups
    # from person.parameter.eRm
      dichX <- ifelse(is.na(X),1,0)
      strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
      gmemb <- as.vector(data.matrix(data.frame(strdata)))
    gindx<-unique(gmemb)
    rsum.all<-rowSums(X,na.rm=T)
    grmed<-tapply(rsum.all,gmemb,median)      #sorted
    ngr<-table(gmemb)                         #sorted
    m.all<-rep(grmed,ngr)                     #sorted,expanded
    rsum.all<-rsum.all[order(gmemb)]
    spl<-ifelse(rsum.all<=m.all,1,2)
    splitcr<-spl
    object$X<-X[order(gmemb),]
  }
}


if (is.numeric(splitcr)){
  if (length(table(splitcr)) > 2) stop("Dichotomous person split required!")
  if (length(splitcr) != dim(object$X)[1]) {
    stop("Mismatch between length of split vector and number of persons!")
  } else {
    rvind <- splitcr
    Xlist <- by(object$X,rvind, function(x) x)
    names(Xlist) <- as.list(sort(unique(splitcr)))
  }}

if (!is.numeric(splitcr)) {
  if (splitcr=="median") {                                   #median split
    rv <- apply(object$X,1,sum,na.rm=TRUE)
    rvsplit <- median(rv)
    rvind <- rep(0,length(rv))
    rvind[rv > rvsplit] <- 1                                 #group with high raw score object
    Xlist <- by(object$X,rvind,function(x) x)
    names(Xlist) <- list("low","high")
    }

  if (splitcr=="mean") {                                     #mean split
    rv <- apply(object$X,1,sum,na.rm=TRUE)
    rvsplit <- mean(rv)
    rvind <- rep(0,length(rv))
    rvind[rv > rvsplit] <- 1                                 #group with highraw scoobject
    Xlist <- by(object$X,rvind,function(x) x)
    names(Xlist) <- list("low","high")
    }

}

del.pos.l <- lapply(Xlist, function(x) {
                    it.sub <- datcheck.LRtest(x,object$X,object$model)  #items to be removed within subgroup
                    })

del.pos <- unique(unlist(del.pos.l))
if ((length(del.pos)) >= (dim(object$X)[2]-1)) {
  stop("\nNo items with appropriate response patterns left to perform Wald-test!\n")
}

if (length(del.pos) > 0) {
    warning("\nThe following items were excluded due to inappropriate response patterns within subgroups: ",immediate.=TRUE)
    cat(colnames(object$X)[del.pos], sep=" ","\n")
    cat("Subgroup models are estimated without these items!\n")
}

if (length(del.pos) > 0) {
  X.el <- object$X[,-(del.pos)]
} else {
  X.el <- object$X
}
Xlist.n <- by(X.el,rvind,function(y) y)
names(Xlist.n) <- names(Xlist)

if (object$model=="RM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg,objectg$betapar,objectg$se.beta)
                               })
       }
if (object$model=="PCM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- PCM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg,objectg$betapar,objectg$se.beta)
                               })
       }
if (object$model=="RSM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RSM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg,objectg$betapar,objectg$se.beta)
                               })
       }

#etapar1 <- unlist(likpar[1,1])
#se1 <- unlist(likpar[2,1])
#etapar2 <- unlist(likpar[1,2])
#se2 <- unlist(likpar[2,2])
#if (length(etapar1) != length(etapar2)) stop("Wald test cannot be performed since number of response item-categories differ over subgroups! \n")

##if (object$model == "RM") {                #for RM-print, which beta are not fixed
##  betapar1 <- likpar[5,1][[1]]
##  betapar2 <- likpar[5,2][[1]]
##  betalab <- colnames(Xlist.n[[1]])        #corresponding item labels
##  hes1 <- likpar[3,1][[1]]                 #Hessian matrices
##  hes2 <- likpar[3,2][[1]]
##  W1 <- likpar[4,1][[1]]                   #design matrixces
##  W2 <- likpar[4,2][[1]]
##  beta1.se <- sqrt(diag(W1%*%solve(hes1)%*%t(W1))) #standard errors for beta1
##  beta2.se <- sqrt(diag(W2%*%solve(hes2)%*%t(W2))) #standard errors for beta2
##} else {
##  betalab <- NULL
##  betapar1 <- NULL
##  betapar2 <- NULL
##  beta1.se <- NULL
##  beta2.se <- NULL
##}

#if (!is.null(betalab)) betalab1 <- betalab[-1]
betapar1 <- likpar[3,][[1]]
beta1.se <- likpar[4,][[1]]
betapar2 <- likpar[3,][[2]]
beta2.se <- likpar[4,][[2]]
num <- (betapar1-betapar2)
denom <- sqrt(beta1.se^2 + beta2.se^2)
W.i <- num/denom
pvalues <- (1-pnorm(abs(W.i)))*2

coef.table <- cbind(W.i,pvalues)
dimnames(coef.table) <- list(names(betapar1),c("z-statistic","p-value"))

result <- list(coef.table=coef.table,betapar1=betapar1,se.beta1=beta1.se,betapar2=betapar2,
se.beta2=beta2.se)
class(result) <- "wald"
result
}

