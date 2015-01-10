\encoding{UTF-8}
\name{RSM}
\alias{RSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of rating scale models}
\description{
  This function computes the parameter estimates of a rating scale model for polytomous
  item responses by using CML estimation. }
\usage{
RSM(X, W, se = TRUE, sum0 = TRUE, etaStart)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix or data frame with item responses (starting from 0); rows represent individuals, columns represent items. Missing values are inserted as \code{NA}.}
  \item{W}{Design matrix for the RSM. If omitted, the function will compute W automatically.}
  \item{se}{If \code{TRUE}, the standard errors are computed.}
  \item{sum0}{If \code{TRUE}, the parameters are normed to sum-0 by specifying
  an appropriate \code{W}. If \code{FALSE}, the first parameter is restricted to 0.}
  \item{etaStart}{A vector of starting values for the eta parameters can be specified. If missing, the 0-vector is used.}
}
\details{
  The design matrix approach transforms the RSM into a partial credit model
  and estimates the corresponding basic parameters by using CML.
  Available methods for RSM-objects are \code{print}, \code{coef}, \code{model.matrix},
  \code{vcov}, \code{summary}, \code{logLik}, \code{person.parameters}, \code{plotICC}, \code{LRtest}.
}
\value{
  Returns an object of class \code{'Rm'}, \code{'eRm'} and contains the log-likelihood value,
  the parameter estimates and their standard errors.

  \item{loglik}{Conditional log-likelihood.}
  \item{iter}{Number of iterations.}
  \item{npar}{Number of parameters.}
  \item{convergence}{See \code{code} output in \code{\link{nlm}}.}
  \item{etapar}{Estimated basic item difficulty parameters (item and category parameters).}
  \item{se.eta}{Standard errors of the estimated basic item parameters.}
  \item{betapar}{Estimated item-category (easiness) parameters.}
  \item{se.beta}{Standard errors of item parameters.}
  \item{hessian}{Hessian matrix if \code{se = TRUE}.}
  \item{W}{Design matrix.}
  \item{X}{Data matrix.}
  \item{X01}{Dichotomized data matrix.}
  \item{call}{The matched call.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations,
Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The \pkg{eRm} package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the \pkg{eRm} package in R. Psychology Science, 49, 26-43.
}
\author{Patrick Mair, Reinhold Hatzinger}
%\note{}

\seealso{\code{\link{RM}},\code{\link{PCM}},\code{\link{LRtest}}
}
\examples{
##RSM with 10 subjects, 3 items
res <- RSM(rsmdat)
res
summary(res)                            #eta and beta parameters with CI
thresholds(res)                         #threshold parameters
}
\keyword{models}
