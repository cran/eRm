`Waldtest.Rm` <-
function(object, splitcr="median")
{
# performs item-based Wald test (Fischer & Molenaar, p.90)
# object... object of class RM
# splitcr... splitting criterion for LR-groups. "median" to a median raw score split,
#            "mean" corobjectponds to the mean raw score split.
#            optionally also a vector of length n for group split can be submitted.


if (is.numeric(splitcr)){
  if (length(table(splitcr)) > 2) stop("Dichotomous person split required!")
  if (length(splitcr) != dim(object$X)[1]) {
    stop("Mismatch between length of split vector and number of persons!")
  } else {
    Xlist <- by(object$X,splitcr, function(x) x)
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

Xlist0 <- lapply(Xlist,function(x) {                                 #eliminate complete NA objectponses
                         tfvec <- apply(x,2,function(z) {
                                               !all(is.na(z))})
                         x[,tfvec]})


#Xlist.n <- lapply(Xlist0,function(x) {       

itdel <- unique(unlist(lapply(Xlist0,function(x) {                         #items which will be deleted
               ri <- apply(x,2,sum,na.rm=TRUE)                       #item raw scoobject
               n.NA <- colSums(apply(x,2,is.na))                     #number of NA's per column
               maxri <- (dim(x)[1]*(apply(x,2,max,na.rm=TRUE)))-n.NA #maximum item raw scoobject with NA
               TFcol <- ((ri==maxri) | (ri==0))                      #full and 0 item scoobject as TRUE
               if (any(TFcol)) {
                 options(warn=0)
                 cat("The following items are not tested since excluded due to 0/full raw score in subgroups:",colnames(x)[TFcol],"\n")
               }
               return(colnames(x)[TFcol])
              })))
              
itkeep <- !(colnames(object$X) %in% itdel)                            #boolean which items should be kept
if (sum(itkeep) == 0) stop("No items left to estimate!")
Xlist.n <- lapply(Xlist0, function(x) x[,itkeep])
       
n.eta <- object$npar

if (object$model=="RM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg)
                               })
       }
if (object$model=="PCM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- PCM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg)
                               })
       }
if (object$model=="RSM") {
       likpar <- sapply(Xlist.n,function(x) {                       #matrix with loglik and npar for each subgroup
                               objectg <- RSM(x)
                               parg <- objectg$etapar
                               seg <- objectg$se.eta
                               list(parg,seg)
                               })
       }

etapar1 <- unlist(likpar[1,1])
se1 <- unlist(likpar[2,1])
etapar2 <- unlist(likpar[1,2])
se2 <- unlist(likpar[2,2])
if (length(etapar1) != length(etapar2)) stop("Wald test cannot be performed since number of response item-categories differ over subgroups! \n")

if (object$model == "RM") {                #for RM-print, which beta are not fixed
  betalab <- colnames(Xlist.n[[1]])        #corresponding item labels
} else {
  betalab <- NULL
}

if (!is.null(betalab)) betalab1 <- betalab[-1]

num <- (etapar1-etapar2)                 #numerator in Wald formula
denom <- sqrt(se1^2 + se2^2)             #denominator
W.i <- (num/denom)                       #z-values
pvalues <- 1-pnorm(abs(W.i))
#pvalues <- 1-pchisq(W.i,1)
   
if (!is.null(betalab)) {                     #for Rasch models item labels are printed out
 coef.table <- as.data.frame(cbind(betalab1,round(W.i,3),round(pvalues,5)))
 dimnames(coef.table) <- list(paste("eta",1:length(etapar1),sep=""),c("Item","z-value","p-value"))
} else {
 coef.table <- cbind(round(W.i,3),round(pvalues,5))
 dimnames(coef.table) <- list(paste("eta",1:length(etapar1),sep=""),c("z-value","p-value"))
}


result <- list(coef.table=coef.table,etapar1=etapar1,se1=se1,etapar2=etapar2,se2=se2,betalab=betalab)
class(result) <- "wald"
result
}

