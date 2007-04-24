`summary.eRm` <-
function(object,labels = TRUE,...)
{

#labels...whether the item parameters should be labelled

cat("\n\n")
cat("Results of",object$model,"fit: \n")
cat("\n")

cat("Log-likelihood:",object$loglik,"\n")
cat("Number of iterations:",object$iter,"\n")
cat("Number of parameters:",object$npar,"\n")
cat("\n")
cat("AIC:",object$IC$AIC,"\n")
cat("BIC:",object$IC$BIC,"\n")
cat("cAIC:",object$IC$cAIC,"\n")
cat("\n")


X <- object$X
X01 <- object$X01
mt_vek <- apply(X,2,max,na.rm=TRUE) 
if (labels) {                 #label beta parameters
  if (object$model=="RM") {                                   #Rasch model
    betanames <- paste("beta",colnames(object$X)) }
  if ((object$model=="RSM") || (object$model=="PCM")) {       #PCM and RSM
    catnames <- sequence(mt_vek)
    itnames <- rep(colnames(object$X),mt_vek)
    betanames <- paste("beta",paste(itnames,catnames,sep="."))
  }
  if ((object$model=="LRSM") || (object$model=="LPCM") || (object$model=="LLTM")) {
    betanames <- paste("beta",colnames(object$X))
  }
} else {                                                      #LLTM, LRSM, LPCM
  betanames <- paste("beta",1:dim(object$W)[1])
}
  
cat("\nItem Parameters (beta):\n")
coeftable <- as.data.frame(cbind(betanames,round(object$betapar,5)))
colnames(coeftable) <- c("Parameter","Estimate")
print(coeftable)
cat("\n")
}

