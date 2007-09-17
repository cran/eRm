`summary.eRm` <-
function(object,...)
{

#labels...whether the item parameters should be labelled

cat("\n")
cat("Results of",object$model,"estimation: \n")
cat("\n")
cat("Call: ", deparse(object$call), "\n")
cat("\n")

cat("Conditional log-likelihood:",object$loglik,"\n")
cat("Number of iterations:",object$iter,"\n")
cat("Number of parameters:",object$npar,"\n")
cat("\n")

X <- object$X
X01 <- object$X01
mt_vek <- apply(X,2,max,na.rm=TRUE)

ci <- confint(object,"eta")
cat("Basic Parameters (eta) with 0.95 CI:\n")
coeftable <- as.data.frame(cbind(round(object$etapar,3),
                           round(object$se.eta,3),round(ci,3)))
colnames(coeftable) <- c("Estimate","Std. Error","lower CI","upper CI")
rownames(coeftable) <- names(object$etapar)
print(coeftable)


ci <- confint(object,"beta")
cat("\nItem Parameters (beta) with 0.95 CI:\n")
coeftable <- as.data.frame(cbind(round(object$betapar,3),
                           round(object$se.beta,3),round(ci,3)))
colnames(coeftable) <- c("Estimate","Std. Error","lower CI","upper CI")
rownames(coeftable) <- names(object$betapar)
print(coeftable)
cat("\n")
}

