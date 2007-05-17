`person.parameter.eRm` <-
function(object)
# estimation of the person parameters with jml
# object of class eRm 
# se... whether standard errors should be computed
# splineInt... whether spline interpolation should be carried out

{

se <- TRUE
splineInt <- TRUE
options(warn=0)
X <- object$X

max.it <- apply(X,2,max,na.rm=TRUE)                               #maximum item raw score without NA
rp <- rowSums(X,na.rm=TRUE)                                       #person raw scores
maxrp <- apply(X,1,function(x.i) {sum(max.it[!is.na(x.i)])})      #maximum item raw score for person i 
TFrow <- ((rp==maxrp) | (rp==0))
pers.exe <- (1:dim(X)[1])[TFrow]                                  #persons excluded from estimation due to 0/full
#pers.notexe <- (1:dim(X)[1])[!TFrow]
X.dummy <- X
if (length(pers.exe) > 0) X <- X[-pers.exe,] 


if (any(is.na(X))) {
  dichX <- ifelse(is.na(X),1,0)
  strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
  gmemb <- as.vector(data.matrix(data.frame(strdata)))
} else {
  gmemb <- rep(1,dim(X)[1])
}



beta <- object$betapar

if (!is.null(object$ngroups))
  if (object$ngroups > 1) stop("Estimation of person parameters for models with group contrasts not possible!")

if (is.null(object$mpoints))  { mpoints <- 1
} else {mpoints <- object$mpoints}
  
mt_vek <- apply(X,2,max,na.rm=TRUE)             #number of categories - 1 for each item
mt_ind <- rep(1:length(mt_vek),mt_vek)          #index i for items
r.pall <- rowSums(X,na.rm=TRUE)                 #person raw scores

X01 <- object$X01
if (length(pers.exe) > 0) X01 <- X01[-pers.exe,]   #if persons excluded due to 0/full response
rowvec <- 1:(dim(X01)[1])

fitlist <- tapply(rowvec,gmemb,function(rind) {         #list with nlm outputs
    
    if (length(rind) > 1) {
       r.i <- colSums(X[rind,],na.rm=TRUE)          #item raw scores
    } else {                                        #if only one person belongs to raw score group
       r.i <- X[rind,]
       r.i[is.na(r.i)] <- 0
    }
    #r.i <- colSums(object$X[rind,],na.rm=TRUE)          #item raw scores
    r.p <- r.pall[rind]                                 #person raw scores for current NA group
    X01g <- X01[rind,]
    X01beta <- rbind(X01g,beta)
    theta <- rep(0,length(r.p))
    
    #==================== ML routines ===================================
    jml.rasch <- function(theta)
    { 
      ksi <- exp(theta)
      denom <- 1/exp(beta)
      lnL <- sum(r.p*theta)-sum(r.i*beta)-sum(log(1+outer(ksi,denom)))
      -lnL
    }
    
    jml <- function(theta) 
    {
        t1t2.list <- tapply(1:(dim(X01beta)[2]),mt_ind, function(xin) {
                                 #print(xin)
                                 xb <- (t(X01beta)[xin,]) #0/1 responses and beta parameters
                                 if (is.vector(xb)) xb <- t(as.matrix(xb))   #if Rasch 
                                 beta.i <- c(0,xb[,dim(xb)[2]])    #item parameter with 0 
                                 x01.i <- (xb[,1:(dim(xb)[2]-1)])    #0/1 matrix for item i without beta
                                 if (is.vector(x01.i)) x01.i <- t(as.matrix(x01.i))
                                 cat0 <- rep(0,dim(x01.i)[2])
                                 cat0[colSums(x01.i)==0] <- 1
                                 x01.i0 <- rbind(cat0,x01.i)       #appending response vector for 0th category
                                 
                                 ind.h <- 0:(length(beta.i)-1)
                                 theta.h <- ind.h %*% t(theta)
                                 term1 <- (theta.h+beta.i)*x01.i0 
                                 t1.i <- sum(colSums(term1))                        #sum over categories and persons
                                                                  
                                 term2 <- exp(theta.h+beta.i)
                                 t2.i <- sum(log(colSums(term2)))                    #sum over categories and persons
                                                              
                                 return(c(t1.i,t2.i))
                               })
      st1st2 <- colSums(matrix(unlist(t1t2.list),ncol=2,byrow=TRUE))                      #sum term1, term2
      lnL <- st1st2[1]-st1st2[2]
      -lnL
    }
    #==================== end ML routines ================================
        
    #==================== call optimizer =================================
    if (object$model == "RM") {
      fit <- nlm(jml.rasch,theta,hessian=se)
    } else {
      fit <- nlm(jml,theta,hessian=se)
    }
    #fit <- optim(theta,jml,method="BFGS",hessian=TRUE)
    
    #=================== end call optimizer ==============================
    loglik <- -fit$minimum
    niter <- fit$iterations
    thetapar <- fit$estimate
    if (se) { 
      se <- sqrt(diag(solve(fit$hessian))) 
    } else {
      se <- NA
      fit$hessian <- NA }

list(loglik=loglik,niter=niter,thetapar=thetapar,se=se,hessian=fit$hessian)
})


loglik <- NULL
niter <- NULL
npar <- NULL
thetapar <- list(NULL) 
se.theta <- list(NULL)
hessian <- list(NULL)
for (i in 1:length(fitlist)) {
  loglik <- c(loglik,fitlist[[i]]$loglik)
  niter <- c(niter,fitlist[[i]]$niter)
  npar <- c(npar,length(fitlist[[i]]$thetapar))
  thetapar[[i]] <- fitlist[[i]]$thetapar
  se.theta[[i]] <- fitlist[[i]]$se
  hessian[[i]] <- fitlist[[i]]$hessian
}

if (splineInt) {                                           #cubic spline interpolation for missing, 0, full raw scores
  x <- rowSums(X,na.rm=TRUE)
  xlist <- split(x,gmemb)
  pred.list <- mapply(function(xx,yy) {
                       y <- tapply(yy,xx, function(xy) {xy[1]})
                       x <- unique(sort(xx))
                       if ((length(x) > 3) || (length(y) > 3)) {        #otherwise splinereg is not admissible
                         fm1 <- interpSpline(x,y)
                         pred.val <- predict(fm1, 0:sum(max.it))
                       } else {
                         warning("Spline interpolation is not performed! Not enough persons in NA subgroups!")
                         NULL
                       }},xlist,thetapar,SIMPLIFY=FALSE)
  X.n <- object$X
  if (any(sapply(pred.list,is.null)))  pred.list <- NULL                           #no spline interpolation applicable
   
} 
                              
result <- list(X=X.n,X01=object$X01,loglik=loglik,npar=npar,iter=niter,
               betapar=object$betapar,thetapar=thetapar,se.theta=se.theta,
               pred.list=pred.list,hessian=hessian,mpoints=mpoints,
               pers.ex=pers.exe,gmemb=gmemb)
class(result) <- "ppar"
result
}
