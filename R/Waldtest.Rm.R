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

Xlist.n <- lapply(Xlist0,function(x) {                               #submatrices without 0 and full item scoobject eliminated
               ri <- apply(x,2,sum,na.rm=TRUE)                       #item raw scoobject
               n.NA <- colSums(apply(x,2,is.na))                     #number of NA's per column
               maxri <- (dim(x)[1]*(apply(x,2,max,na.rm=TRUE)))-n.NA #maximum item raw scoobject with NA
               TFcol <- ((ri==maxri) | (ri==0))                      #full and 0 item scoobject as TRUE
               if (any(TFcol)) {
                 options(warn=0)
                 
                 cat("The following items are deleted due to 0/full raw score:",colnames(x)[TFcol],"\n")
                 warning("Waldtest may not work, choose another split!") 
               }
               x.n <- x[,TFcol==FALSE]
               if (dim(x.n)[2]==0) x.n <- NULL                        #nothing left to estimate
               x.n
              })

Xlist.n <- Xlist.n[!sapply(Xlist.n,is.null)]
Xmax <- apply(object$X,2,max,na.rm=TRUE)
lapply(Xlist.n,function(x) {
                 submax <- apply(x,2,max,na.rm=TRUE)
                 if (length(submax) != length(Xmax)) {
                   stop("Number of parameters is different for subgroups. Waldtest cannot be performed, choose another split!")
                 } else {  
                   if (any(submax != Xmax)) {
                     stop("Maximal category responses in subgroups do not match. Waldtest cannot be performed!") 
                 }}})

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

if (object$model == "RM") {                #for RM-print, which beta are not fixed
  betalab <- colnames(object$X)[-1] #corresponding item labels
} else {
  betalab <- NULL
}


result <- list(etapar1=etapar1,se1=se1,etapar2=etapar2,se2=se2,betalab=betalab)
class(result) <- "wald"
result
}

