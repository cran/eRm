\encoding{UTF-8}
\name{RM}
\alias{RM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of Rasch Models}
\description{
  This function computes the parameter estimates of a Rasch model for binary item responses by using CML estimation.
}
\usage{
RM(X, W, se = TRUE, sum0 = TRUE, etaStart)
}
\arguments{
  \item{X}{Input 0/1 data matrix or data frame; rows represent individuals, columns represent items. Missing values are inserted as \code{NA}.}
  \item{W}{Design matrix for the Rasch model. If omitted, the function will compute W automatically.}
  \item{se}{If \code{TRUE}, the standard errors are computed.}
  \item{sum0}{If \code{TRUE}, the parameters are normed to sum-0 by specifying
  an appropriate \code{W}. If \code{FALSE}, the first parameter is restricted to 0.}
  \item{etaStart}{A vector of starting values for the eta parameters can be specified. If missing, the 0-vector is used.}
}
\details{
  For estimating the item parameters the CML method is used.
  Available methods for RM-objects are:\cr
  \code{print}, \code{coef}, \code{model.matrix},
  \code{vcov}, \code{summary}, \code{logLik}, \code{person.parameter}, \code{LRtest},
  \code{Waldtest}, \code{plotICC}, \code{plotjointICC}.
}
\value{
  Returns an object of class \code{dRm, Rm, eRm} and contains the log-likelihood value, the parameter estimates and their standard errors.

  \item{loglik}{Conditional log-likelihood.}
  \item{iter}{Number of iterations.}
  \item{npar}{Number of parameters.}
  \item{convergence}{See \code{code} output in \code{\link{nlm}}.}
  \item{etapar}{Estimated basic item difficulty parameters.}
  \item{se.eta}{Standard errors of the estimated basic item parameters.}
  \item{betapar}{Estimated item (easiness) parameters.}
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
\seealso{\code{\link{RSM}},\code{\link{PCM}}, \code{\link{LRtest}}, \code{\link{Waldtest}}
}
\examples{
# Rasch model with beta.1 restricted to 0
res <- RM(raschdat1, sum0 = FALSE)
res
summary(res)
res$W                                       #generated design matrix

# Rasch model with sum-0 beta restriction; no standard errors computed
res <- RM(raschdat1, se = FALSE, sum0 = TRUE)
res
summary(res)
res$W                                       #generated design matrix

#Rasch model with missing values
res <- RM(raschdat2)
res
summary(res)
}
\keyword{models}
