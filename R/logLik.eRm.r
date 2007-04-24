logLik.eRm <- function(object,...)
{
#object of class eRm
  REML <- TRUE
  if (REML)
        stop("cannot calculate REML log-likelihood for \"eRm\" objects")

  val <- object$loglik
  attr(val, "df") <- object$npar
  class(val) <- "logLik"
  val
}


