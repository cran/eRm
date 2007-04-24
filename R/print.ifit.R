`print.ifit` <-
function(x,...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{

  pvalues <- mapply(function(inf,idf) {                                  #p-values for item infit and outfit statistics
                      p.infit <- 1-pchisq(inf,idf)
                      list(p.infit)
                    },x$i.fit,x$i.df,SIMPLIFY=FALSE)

  cat("\nChi-square Itemfit Statistics: \n")
  for (i in 1:length(x$i.fit)) {
    if (length(x$i.fit) > 1) {cat("Person Group:",i,"\n")}
    coef.table <- cbind(x$i.fit[[i]],x$i.df[[i]],pvalues[[i]][[1]])
    dimnames(coef.table) <- list(names(x$i.fit[[1]]),c("Itemfit","df","p.value"))
    print(coef.table)
    cat("\n")
  }
  
}

