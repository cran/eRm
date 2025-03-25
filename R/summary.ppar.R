`summary.ppar` <-
function(object,...)
# summary method for object of class "ppar"
{
  
  if (length(object$pers.ex) > 0) {
    thetaind <- rownames(object$X)[-object$pers.ex]
  } else {
    thetaind <- rownames(object$X)
  }
  
  if (length(object$pers.ex) > 0) {    
    X <- object$X[-object$pers.ex,]                                        #list with raw scores
    sumlist <- by(object$X[-object$pers.ex,],object$gmemb,rowSums,na.rm=TRUE)
  } else {
    X <- object$X
    sumlist <- by(object$X,object$gmemb,rowSums,na.rm=TRUE)
  }
  
    
  if (any(is.na(object$X))) {                                       #recompute gmemb without persons excluded
    dichX <- ifelse(is.na(object$X),1,0)
    strdata <- apply(dichX,1,function(x) {paste(x,collapse="")})
    gmemb <- as.vector(data.matrix(data.frame(strdata)))
  } else {
    gmemb <- rep(1,dim(object$X)[1])
  }
  
  cat("\n")
  cat("Estimation of Ability Parameters")
  for (i in 1:length(object$thetapar)) {
    cat("\n\n")
    if (length(object$thetapar) > 1) {
      cat("Person NA Group:",i,"\n")
      xvec <- rep(NA, (dim(X)[2]))
      notNApos <- which(!is.na(as.vector(rbind(X[object$gmemb == i,])[1,])))
      xvec[notNApos] <- "x"
      cat("NA pattern:",xvec,"\n")
      }
    cat("Collapsed log-likelihood:",object$loglik[[i]],"\n")
    cat("Number of iterations:",object$iter[[i]],"\n")
    cat("Number of parameters:",object$npar[[i]],"\n")
    cat("\n")
    cat("ML estimated ability parameters (without spline interpolated values): \n")
    coef.table <- cbind(object$thetapar[[i]],object$se.theta[[i]],confint(object)[[i]])
    dimnames(coef.table) <- list(paste("theta",thetaind[object$gmemb==i]),c("Estimate","Std. Err.",colnames(confint(object)[[i]])))
    print(coef.table)
  }
}

