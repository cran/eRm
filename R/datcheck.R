"datcheck" <-
function(X,W)
{
  if (is.data.frame(X)==TRUE) X <- as.matrix(X)                  #X as data frame allowed
  if (sum(is.na(X))!= 0) stop("No missing values allowed in X!") #missing value check
  #if (W != NA) {                                                 #check whether W is specified correctly
  #  WX <- (dim(W)[1]==dim(X)[2])
  #  if (WX == FALSE) stop("Wrong specification of W!")
  #}
list(X=X)
}

