`print.eRm` <-
function(x,...)  {                                         #print method for all models
  cat("\n")
  cat("Results of", x$model, "estimation: \n")
  cat("\n")
  cat("Call: ", deparse(x$call), "\n")
  cat("\n")
  cat("Conditional log-likelihood:", x$loglik, "\n")
  cat("Number of iterations:", x$iter, "\n")
  cat("Number of parameters:", x$npar, "\n")
  cat("\n")
  cat("Basic Parameters eta:")                             #eta parameters
  cat("\n")
  etapar <- x$etapar
  #nameeta <- paste("eta",1:dim(x$W)[2])
  se <- x$se.eta
  result <- rbind(etapar, se)
  #colnames(result) <- nameeta
  rownames(result) <- c("Estimate", "Std.Err")
  print(result)
  cat("\n\n")
}

