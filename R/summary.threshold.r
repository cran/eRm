summary.threshold <- function(object,...)
{
#object of class "threshold"

  coef.table <- cbind(object$thresh,object$se.thresh,confint(object))
  dimnames(coef.table) <- list(names(object$thresh),c("Estimate","Std. Err.",colnames(confint(object))))
  cat("\n")
  print(coef.table)
  cat("\n")
}