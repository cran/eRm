`summary.eRm` <-
function(object,...)
{

cat("\n\n")
cat("Results of",object$model,"fit\n")
cat("\n")

cat("number of iterations: ",object$iter,"\n")
cat("log likelihood: ",object$loglik,"\n")
cat("df = ",object$df,"\n")

etanames <- paste("eta",1:length(object$etapar))
#etanames <- names(as.data.frame(object$W))  

cat("\nItem Parameters:\n\n")
tvalue <- object$etapar/object$se_eta
pvalue <- 2 * pnorm(-abs(tvalue))
dn <- c("Estimate", "Std. Error")
coef.table <- cbind(object$etapar, object$se_eta, tvalue, pvalue)
dn <- c("Estimate", "Std. Error")
dimnames(coef.table) <- list(etanames, c(dn,"z value", "Pr(>|z|)"))
print(coef.table)
cat("\n\n")

if (length(object$LR)>1) {
  cat("Goodness-of-fit (modified Andersen's test):\n\n")
  cat("LR statistic =",object$LR$Chisq,"  df =",object$LR$df, " p =",object$LR$pvalue,"\n\n")
}
}

