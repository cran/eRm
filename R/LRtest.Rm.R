`LRtest.Rm` <-
function(object, splitcr="median",alpha=0.05,se=TRUE)
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
    rv <- apply(object$X,1,sum,na.rm=TRUE)                      #person raw scoobject
    Xlist <- by(object$X,rv,function(x) x)
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
    Xlist <- by(object$X,splitcr, function(x) x) 
    names(Xlist) <- as.list(sort(unique(splitcr)))
    }}
    
Xlist0 <- lapply(Xlist,function(x) {                                 #eliminate complete NA objectponses
                         tfvec <- apply(x,2,function(z) {
                                               !all(is.na(z))})
                         x[,tfvec]})

Xlist.n <- lapply(Xlist0,function(x) {                               #submatrices without 0 and full item scoobject eliminated
               x <- as.matrix(x)
               ri <- apply(x,2,sum,na.rm=TRUE)                       #item raw scoobject
               n.NA <- colSums(as.matrix((apply(x,2,is.na))))        #number of NA's per column
               maxri <- (dim(x)[1]*(apply(x,2,max,na.rm=TRUE)))-n.NA #maximum item raw scoobject with NA
               TFcol <- ((ri==maxri) | (ri==0))                      #full and 0 item scoobject as TRUE
               x.n <- x[,TFcol==FALSE]
               if (length(x.n) == 0) {
                 x.n <- NULL 
               } else {
                 if (dim(x.n)[2]==0) x.n <- NULL                     #nothing left to estimate
               }
               x.n
              })

Xlist.n <- Xlist.n[!sapply(Xlist.n,is.null)]                        #delete NULL elements in Xlist.n


if (object$model=="RM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               betalab <- colnames(objectg$X)
                               list(likg,nparg,objectg$betapar,betalab,objectg$etapar,objectg$se.eta)
                               })
       }
if (object$model=="PCM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- PCM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,NULL,objectg$etapar,objectg$se.eta)
                               })
       }
if (object$model=="RSM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RSM(x,se=se)
                               likg <- objectg$loglik
                               nparg <- length(objectg$etapar)
                               list(likg,nparg,objectg$betapar,NULL,objectg$etapar,objectg$se.eta)
                               })
       }
       
loglikg <- sum(unlist(likpar[1,]))                    #sum of likelihood value for subgroups
LR <- 2*(abs(loglikg-object$loglik))                  #LR value
df = sum(unlist(likpar[2,]))-(length(object$etapar))  #final degrees of freedom
Chisq <- qchisq(1-alpha,df)
pvalue <- 1-pchisq(LR,df)                             #pvalue

betalist <- likpar[3,]                                #organizing betalist
#betalist <- lapply(lapply(betalist,as.matrix),t)
if (object$model == "RM") {                           #label betalist
  betalab <- likpar[4,] 
  for (i in 1:length(betalist)) names(betalist[[i]]) <- betalab[[i]] 
} else {
  for (i in 1:length(betalist)) names(betalist[[i]]) <- 1:length(betalist[[i]])
} 
         
result <- list(X=object$X, X.list=Xlist.n, model=object$model,LR=LR,
               Chisq=Chisq, df=df, pvalue=pvalue, likgroup=unlist(likpar[1,],use.names=FALSE),
               betalist=betalist, etalist=likpar[5,],selist=likpar[6,])
class(result) <- "LR"
result
}

