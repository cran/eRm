\encoding{UTF-8}
\name{print.eRm}
\alias{print.eRm}
\alias{summary.eRm}
\alias{vcov.eRm}
\alias{model.matrix.eRm}
\alias{coef.eRm}
\alias{logLik.eRm}
\alias{confint.eRm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for extended Rasch models}
\description{Several methods for objects of class \code{'eRm'}.}
\usage{
\method{print}{eRm}(x, ...)
\method{summary}{eRm}(object, ...)
\method{coef}{eRm}(object, parm="beta", ...)
\method{model.matrix}{eRm}(object, ...)
\method{vcov}{eRm}(object, ...)
\method{logLik}{eRm}(object, ...)
\method{confint}{eRm}(object, parm = "beta", level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{eRm}.}
  \item{object}{Object of class \code{eRm}.}
  \item{parm}{Either \code{"eta"} or \code{"beta"}.}
  \item{level}{Alpha-level.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}
}
\details{
  The \code{print} method displays  the value of
  the log-likelihood, parameter estimates (basic parameters eta) and their standard errors.
  For RM, RSM, and PCM models, the etas are difficulty parameters, for the LLTM, LRSM,
  LPCM the sign of the parameters depend on the design matrix and are easiness effects by default.
  The \code{summary} method additionally gives the full set of item parameters beta as
  easiness parameters for all models.

  Print methods are also available for the functions \code{logLik} and \code{confint}
  (see below).
}

\value{
The methods below are extractor functions and return various quantities:
  \code{vcov} returns the variance-covariance matrix of the parameter estimates,
  \code{coef} a vector of estimates of the eta or beta basic parameters,
  \code{model.matrix} the design matrix,
  \code{logLik} an object with elements \code{loglik} and \code{df} containing
  the log-likelihood value and df.
  \code{confint} a matrix of confidence interval values for eta or beta.
}
\author{Patrick Mair, Reinhold Hatzinger}

\examples{
res <- RM(raschdat1)
res
summary(res)
coef(res)
vcov(res)
model.matrix(res)
logLik(res)
}
\keyword{models}

