\name{print.eRm}
\alias{print.eRm}
\alias{summary.eRm}
\alias{vcov.eRm}
\alias{model.matrix.eRm}
\alias{coef.eRm}
\alias{logLik.eRm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for extended Rasch models}
\description{}
\usage{
\method{print}{eRm}(x, ...)
\method{summary}{eRm}(object, labels = TRUE, ...)
\method{coef}{eRm}(object, ...)
\method{model.matrix}{eRm}(object, ...)
\method{vcov}{eRm}(object, ...)
\method{logLik}{eRm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{eRm}.}
  \item{object}{Object of class \code{eRm}.}
  \item{labels}{If \code{TRUE}, items are labeled.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}
}
\details{
  \code{vcov} returns the variance-covariance matrix of the parameter estimates, 
  \code{coef}   the estimates of the basic parameters, \code{print} the value of 
  the log-likelihood and the   parameter estimates and their standard errors, 
  \code{model.matrix} the design matrix, \code{logLik} the log-likelihood value.
}
\author{Patrick Mair, Reinhold Hatzinger}

\examples{
data(raschdat1)
res <- RM(raschdat1)
print(res)
summary(res)
coef(res)
vcov(res)
model.matrix(res)
}

\keyword{models}

