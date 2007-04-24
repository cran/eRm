`print.wald` <-
function(x,...)
#print method for objects of class "wald" (from waldtest)
{

   num <- (x$etapar1-x$etapar2)                 #numerator in Wald formula
   denom <- sqrt(x$se1^2 + x$se2^2)                     #denominator
   W.i <- (num/denom)                         #z-values
   pvalues <- 1-pnorm(W.i)
   #pvalues <- 1-pchisq(W.i,1)
   
   if (!is.null(x$betalab)) {                     #for Rasch models item labels are printed out
     coef.table <- as.data.frame(cbind(x$betalab,round(W.i,3),round(pvalues,5)))
     dimnames(coef.table) <- list(paste("eta",1:length(x$etapar1),sep=""),c("item","z-value","p-value"))
   } else {
     coef.table <- cbind(round(W.i,3),round(pvalues,5))
     dimnames(coef.table) <- list(paste("eta",1:length(x$etapar1),sep=""),c("z-value","p-value"))
   }
   
   cat("\nWald test on item level (z-values):\n\n")
   print(coef.table)
   cat("\n")
}

