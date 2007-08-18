`print.pfit` <-
function(x,...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{
  pvalues <- 1-pchisq(x$p.fit,x$p.df)
  coef.table <- cbind(round(x$p.fit,3),x$p.df,round(pvalues,3))
  colnames(coef.table) <- c("Chisq","df","p-value")
  rownames(coef.table) <- names(x$p.fit)
  cat("\nChi-square Personfit Statistics: \n")
  print(coef.table)
  cat("\n")
}

