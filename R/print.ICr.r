print.ICr <- function(x,...)
{
#print method for objects of class "ICr" (from function "IC")

cat("\nUnconditional (joint) log-likelihood:",x$j.loglik)
cat("\n")
cat("AIC:",x$AIC,"\n")
cat("BIC:",x$BIC,"\n")
cat("cAIC:",x$cAIC,"\n\n")
}