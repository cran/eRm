`print.eRm` <-
function(x,...)  {                                         #print method for all models
  cat("\n")
  ll <- x$loglik
  cat("log-likelihood: ",ll)
  cat("\n")
  cat("\n")
  etapar <- x$etapar
  
  nameeta <- paste("eta",1:dim(x$W)[2])
  #nameeta <- names(as.data.frame(x$W))                    
  
  se <- x$se_eta
  result <- rbind(etapar, se)                                         #(virtual) item parameters, standard errors
  
  colnames(result) <- nameeta
  rownames(result) <- c("Estimate", "Std.Err")
  print(result)
  cat("\n")
  invisible(x)
}

