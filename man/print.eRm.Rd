\name{print.eRm}
\alias{plot.Rm}
\alias{print.eRm}
\alias{summary.Rm}
\alias{vcov.eRm}
\alias{model.matrix.eRm}
\alias{coef.eRm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods for extended Rasch models}
\description{The \code{summary} and \code{plot} methods produce LR-tests and goodness-of-fit
plots for objects of class \code{Rm}, i.e. the Rasch model, the RSM, and the PCM. The remaining methods
can be used for any type of eRm-models.
}
\usage{
\method{print}{eRm}(x, ...)
\method{plot}{Rm}(x,y, ...)
\method{summary}{Rm}(object, ...)
\method{coef}{eRm}(object, ...)
\method{model.matrix}{eRm}(object, ...)
\method{vcov}{eRm}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x, y}{An object of class \code{eRm} or \code{Rm}.}
  \item{\dots}{arguments to pass down (ignored).}
}
\details{
  \code{vcov} returns the variance-covariance matrix of the parameter estimates, \code{coef} 
  the estimates of the basic parameters, \code{print} the value of the log-likelihood and the
  parameter estimates and their standard errors, \code{model.matrix} the design matrix, and 
  \code{summary} and \code{plot} as described above.
}
\author{Patrick Mair, Reinhold Hatzinger}

\examples{
data(X_rasch)
res <- RM(X_rasch)
print(res)
plot(res)
summary(res)
coef(res)
vcov(res)
model.matrix(res)
}

\keyword{models}

