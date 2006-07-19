"lr" <-
function(likall)
{
  likvek <- sapply(likall,function(lg) {lik <- (lg[[1]]$minimum)}) #get out vektor of likelihood values
                              
  LR <- abs(2*(sum(likvek[2:length(likvek)])-likvek[1]))           #LR statistic - Andersen (1972)
  n_eta <- length(likall[[1]][[1]]$estimate)                       #number of eta parameters
  df = 2*n_eta-n_eta                                          #2 because of two subgroups
  pvalue <- 1-pchisq(LR,df)
  list(Chisq=LR,df=df,pvalue=pvalue)
}

