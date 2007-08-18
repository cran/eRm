`print.eRm` <-
function(x,...)  {                                         #print method for all models
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
                                                         #(virtual) item parameters beta
  cat("Item Parameters beta:")
  cat("\n")
  betapar <- x$betapar
  se <- x$se.beta
  result <- rbind(betapar, se)                                         
  rownames(result) <- c("Estimate", "Std.Err.")
  print(result)
  cat("\n\n")
}

