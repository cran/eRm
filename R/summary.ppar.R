`summary.ppar` <-
function(object,...)
# summary method for object of class "ppar"
{
  
  if (length(object$pers.ex) > 0) {
    thetaind <- rownames(object$X)[-object$pers.ex]
  } else {
    thetaind <- rownames(object$X)
  }
    
  cat("\n")
  cat("Estimation of Person Parameters")
  for (i in 1:length(object$thetapar)) {
    cat("\n\n")
    if (length(object$thetapar) > 1) {cat("Person Group:",i,"\n")}
    cat("Log-likelihood:",object$loglik[[i]],"\n")
    cat("Number of iterations:",object$iter[[i]],"\n")
    cat("Number of parameters:",object$npar[[i]],"\n")
    cat("\n")
    cat("ML estimated person parameters (without spline approximations): \n")
    coef.table <- cbind(object$thetapar[[i]],object$se.theta[[i]])
    dimnames(coef.table) <- list(paste("theta",thetaind[object$gmemb==i]),c("Estimate","Std. Error."))
    print(coef.table)
  }
}

