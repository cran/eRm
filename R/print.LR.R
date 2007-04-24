`print.LR` <-
function(x,...)
{
#print method for object of class "LR" (LRtest)
  cat("\n")
  cat("Andersen LR-test: \n")
  cat("LR-value:", x$LR,"\n")
  cat("df:",x$df,"\n")
  cat("p-value: ",x$pvalue,"\n")
  cat("\n")
}

