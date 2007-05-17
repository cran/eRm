`print.pfit` <-
function(x,...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{

  pvalues <- mapply(function(inf,idf) {                                  #p-values for item infit and outfit statistics
                      p.infit <- 1-pchisq(inf,idf)
                      list(p.infit)
                    },x$p.fit,x$p.df,SIMPLIFY=FALSE)

  cat("\nChi-square Personfit Statistics: \n")
  for (i in 1:length(x$p.fit)) {
    if (length(x$p.fit) > 1) {cat("Person NA Group:",i,"\n")}
    coef.table <- cbind(x$p.fit[[i]],x$p.df[[i]],pvalues[[i]][[1]])
    dimnames(coef.table) <- list(names(x$p.fit[[i]]),c("Personfit","df","p.value"))
    print(coef.table)
    cat("\n")
  }
  
}

