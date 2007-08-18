`LRtest.Rm` <-
function(object, splitcr="median", se=FALSE)
{
# performs Andersen LR-test
# object... object of class RM
# splitcr... splitting criterion for LR-groups. "all.r" corresponds to a complete
#            raw score split (r=1,...,k-1), "median" to a median raw score split,
#            "mean" corresponds to the mean raw score split. 
#            optionally also a vector of length n for group split can be submitted.
# se...whether standard errors should be computed

if (!is.numeric(splitcr)) {
  if (splitcr=="all.r") {                                    #full raw score split
    rvind <- apply(object$X,1,sum,na.rm=TRUE)                      #person raw scoobject
    Xlist <- by(object$X,rvind,function(x) x)
    names(Xlist) <- as.list(sort(unique(rv)))
    }

  if (splitcr=="median") {                                   #median split
    cat("Warning message: Persons with median raw scores are assigned to the lower raw score group!\n")
    rv <- apply(object$X,1,sum,na.rm=TRUE)
    rvsplit <- median(rv)
    rvind <- rep(0,length(rv))
    rvind[rv > rvsplit] <- 1                                 #group with highraw scoobject
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

if (is.numeric(splitcr)) {                                 #manual raw score split
  if (length(splitcr)!=dim(object$X)[1]){
    stop("Mismatch between length of split vector and number of persons!")
  } else {
    rvind <- splitcr
    Xlist <- by(object$X,rvind, function(x) x) 
    names(Xlist) <- as.list(sort(unique(splitcr)))
    }}
    
#----------item to be deleted---------------
del.pos.l <- lapply(Xlist, function(x) {
                    it.sub <- datcheck.LRtest(x,object$X,object$model)  #items to be removed within subgroup
                    })

del.pos <- unique(unlist(del.pos.l)) 
if ((length(del.pos)) >= (dim(object$X)[2]-1)) {
  stop("\nNo items with appropriate response patterns left to perform LR-test!\n")
}

if (length(del.pos) > 0) {
  warning("\nThe following items were excluded due to inappropriate response patterns within subgroups: ",immediate.=TRUE) 
    cat(colnames(object$X)[del.pos], sep=" ","\n")
    cat("Full and subgroup models are estimated without these items!\n")
}
            

if (length(del.pos) > 0) {
  X.el <- object$X[,-(del.pos)]
} else {
  X.el <- object$X
}                    
Xlist.n <- by(X.el,rvind,function(y) y)
names(Xlist.n) <- names(Xlist)
if (length(del.pos) > 0) Xlist.n <- c(Xlist.n,list(X.el))
              
if (object$model=="RM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                              # betalab <- colnames(objectg$X)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)
                               })
       }
if (object$model=="PCM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- PCM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)
                               })
       }
if (object$model=="RSM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RSM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,objectg$etapar,objectg$se.beta)
                               })
       }
       
if (length(del.pos) > 0) {                  #re-estimate full model
  pos <- length(Xlist.n)                    #position of the full model
  loglik.all <- likpar[1,pos][[1]]          #loglik full model
  etapar.all <- rep(0,likpar[2,pos])         #etapar full model (filled with 0 for df computation)
  likpar <- likpar[,-pos]
  Xlist.n <- Xlist.n[-pos]
} else {
  loglik.all <- object$loglik
  etapar.all <- object$etapar 
}     
       
loglikg <- sum(unlist(likpar[1,]))                    #sum of likelihood value for subgroups
LR <- 2*(abs(loglikg-loglik.all))                  #LR value
df = sum(unlist(likpar[2,]))-(length(etapar.all))  #final degrees of freedom
pvalue <- 1-pchisq(LR,df)                             #pvalue

betalist <- likpar[3,]                                #organizing betalist
         
result <- list(X=object$X, X.list=Xlist.n, model=object$model,LR=LR,
               df=df, pvalue=pvalue, likgroup=unlist(likpar[1,],use.names=FALSE),
               betalist=betalist, etalist=likpar[4,],selist=likpar[5,])
class(result) <- "LR"
result
}

