`print.ifit` <-
function(x, visible = TRUE, sort_by = c("none", "p", "outfit_MSQ", "infit_MSQ", "outfit_t", "infit_t", "discrim"), 
         decreasing = FALSE, digits = 3, ...)
# print method for itemfit
# x...object of class "ifit" from (itemfit)
{
  sort_by <- match.arg(sort_by,  c("none", "p", "outfit_MSQ", "infit_MSQ", "outfit_t", "infit_t", "discrim"), several.ok = FALSE)
  
  pvalues <- 1-pchisq(x$i.fit,x$i.df-1)  
  coef.table <- cbind(round(x$i.fit,digits), x$i.df-1, round(pvalues,digits), round(x$i.outfitMSQ, digits), 
                      round(x$i.infitMSQ, digits),round(x$i.outfitZ, digits),round(x$i.infitZ, digits),
                      round(x$i.disc, digits))
  colnames(coef.table) <- c("Chisq","df","p-value","Outfit MSQ", "Infit MSQ", "Outfit t", "Infit t", "Discrim")
  rownames(coef.table) <- names(x$i.fit)
  
  ## sort 
  if (length(sort_by) == 6) sort_by <- "none"
  if (sort_by != "none") {
    sort_by <- sort_by[1]
    if (sort_by == "p") sortcol <- 3
    if (sort_by == "outfit_MSQ") sortcol <- 4
    if (sort_by == "infit_MSQ") sortcol <- 5
    if (sort_by == "outfit_t") sortcol <- 6
    if (sort_by == "infit_t") sortcol <- 7
    if (sort_by == "discrim") sortcol <- 8
    if (sort_by != "discrim") {
      ind <- order(abs(coef.table[, sortcol]), decreasing = decreasing)
    } else {
      ind <- order(coef.table[, sortcol], decreasing = decreasing)
    }  
    coef.table <- coef.table[ind,]
  }
  
  if (visible){       
    cat("\nItemfit Statistics: \n")
    print(coef.table)
    cat("\n")
  }
  invisible(coef.table)
}

