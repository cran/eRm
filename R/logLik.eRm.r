logLik.eRm <- function(object,...)
{
#object of class eRm
  val <- object$loglik
  attr(val, "df") <- object$npar
  class(val) <- "logLik.eRm"
  val
}


