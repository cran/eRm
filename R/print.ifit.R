`print.ifit` <-
function(x,...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{
  pvalues <- 1-pchisq(x$i.fit,x$i.df)
  coef.table <- cbind(round(x$i.fit,3),x$i.df,round(pvalues,3))
  colnames(coef.table) <- c("Chisq","df","p-value")
  rownames(coef.table) <- names(x$i.fit)
  cat("\nChi-square Itemfit Statistics: \n")
  print(coef.table)
  cat("\n")
}

