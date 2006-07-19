"print.eRm" <-
function(x,...)  {                                         #print method for all models
  cat("\n")
  ll <- x$loglik
  cat("log-likelihood: ",ll)
  cat("\n")
  cat("\n")
  etapar <- coef(x)
  etapar1 <- c(0,etapar)                                                #include parameter item 1 set to 0
  Tcontr <- x$mpoints-1                                                 #number of main and interaction effects
  Gcontr <- x$ngroups-1
  TGcontr <- Tcontr*Gcontr
  nitems <- length(etapar1)-Tcontr-Gcontr-TGcontr                       #number of item parameters
  
  se <- x$se_eta
  se1 <- c(NA,se)
  result <- rbind(etapar1, se1)                                         #(virtual) item parameters, standard errors
  
  item_n <- paste("Item",1:nitems,sep="")
  if (x$mpoints == 1) {nmT <- NULL
  }else nmT <- paste("T",2:x$mpoints,sep="")
  if (x$ngroups == 1) {nmG <- NULL
  }else nmG <- paste("G",2:x$ngroups,sep="")
  if ((x$mpoints > 1) && (x$ngroups > 1)) {nwTG <- paste(nmT,":",nmG,sep="") 
  }else nwTG <- NULL 
  
  tg_str <- c(nmT,nmG,nwTG)                                       #main and interaction effects
  all_str <- c(item_n,tg_str)
  if (length(all_str) == dim(result)[2]) {colnames(result) <- all_str
  }else {
     tg_str <- paste("e",1:(length(etapar1)-nitems),sep="")
     colnames(result) <- c(item_n,tg_str)
  }
  
  rownames(result) <- c("Estimate", "Std.Err")
  print(result)
  cat("\n")
  invisible(x)
}

