logLik.ppar <- function(object,...)
{
#object of class ppar
  val <- IC(object)$j.loglik
  attr(val, "df") <- dim(object$W)[2]
  class(val) <- "logLik.ppar"
  val
}