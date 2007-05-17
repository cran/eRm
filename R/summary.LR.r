summary.LR <- function(object,...)
# summary method for objects of class "LR" (from LRtest")
{
  cat("\n")
  cat("Andersen LR-test: \n")
  cat("LR-value:", object$LR,"\n")
  cat("Critical Chi-squared value: ",object$Chisq,"\n")
  cat("df:",object$df,"\n")
  cat("p-value: ",object$pvalue,"\n")
  cat("\n")
  
  mt_vek <- apply(object$X,2,max,na.rm=TRUE) 
  #if ((object$model=="RSM") || (object$model=="PCM")) {       #PCM and RSM
  #  catnames <- sequence(mt_vek)
  #  itnames <- rep(colnames(object$X),mt_vek)
  #  betanames <- paste("beta",paste(itnames,catnames,sep="."))
  #}
  
  
  for (i in 1:length(object$betalist)) {
    cat("\n")
    cat("Splitted person subgroup",i)
    cat("\n")
    cat("Log-likelihood: ",object$likgroup[i])
    cat("\n\n")
    cat("Beta Parameters: \n")
    betavec <- object$betalist[[i]]
    #if (is.null(names(betavec))) names(betavec) <- betanames
    print(betavec)
    cat("\n")
  }
}

    