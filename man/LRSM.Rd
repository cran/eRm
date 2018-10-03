\encoding{UTF-8}
\name{LRSM}
\alias{LRSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of linear rating scale models}
\description{
  This function computes the parameter estimates of a linear rating scale model (LRSM)
  for polytomuous item responses by using CML estimation.
}
\usage{
LRSM(X, W , mpoints = 1, groupvec = 1, se = TRUE, sum0 = TRUE,
   etaStart)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix or data frame; rows represent individuals (N in total), columns represent items. Missing values are inserted as \code{NA}.}
  \item{W}{Design matrix for the LRSM. If omitted, the function will compute W automatically.}
  \item{mpoints}{Number of measurement points.}
  \item{groupvec}{Vector of length N which determines the group membership of each subject, starting from 1}
  \item{se}{If \code{TRUE}, the standard errors are computed.}
  \item{sum0}{If \code{TRUE}, the parameters are normalized to sum-0 by specifying
  an appropriate \code{W}. If \code{FALSE}, the first parameter is restricted to 0.}
  \item{etaStart}{A vector of starting values for the eta parameters can be specified. If missing, the 0-vector is used.}
}
\details{
  Through appropriate definition of \code{W} the LRSM can be viewed as a more parsimonous
  RSM, on the one hand, e.g. by imposing some cognitive base operations
  to solve the items. One the other hand, linear extensions of the Rasch model
  such as group comparisons and repeated measurement designs can be computed.
  If more than one measurement point is examined, the item responses for the 2nd, 3rd, etc.
  measurement point are added column-wise in X.

  If \code{W} is user-defined, it is nevertheless necessary to
  specify \code{mpoints} and \code{groupvec}. It is important that first the time contrasts and
  then the group contrasts have to be imposed.

  Available methods for LRSM-objects are:
  \code{print}, \code{coef},
  \code{model.matrix}, \code{vcov},\code{summary}, \code{logLik}, \code{person.parameters}.
}
\value{
  Returns an object of class \code{'eRm'} containing:

  \item{loglik}{Conditional log-likelihood.}
  \item{iter}{Number of iterations.}
  \item{npar}{Number of parameters.}
  \item{convergence}{See \code{code} output in \code{\link{nlm}}.}
  \item{etapar}{Estimated basic item parameters (item and category parameters).}
  \item{se.eta}{Standard errors of the estimated basic item parameters.}
  \item{betapar}{Estimated item (easiness) parameters.}
  \item{se.beta}{Standard errors of item parameters.}
  \item{hessian}{Hessian matrix if \code{se = TRUE}.}
  \item{W}{Design matrix.}
  \item{X}{Data matrix.}
  \item{X01}{Dichotomized data matrix.}
  \item{groupvec}{Group membership vector.}
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
\seealso{\code{\link{LLTM}},\code{\link{LPCM}}}
\examples{
#LRSM for two measurement points
#20 subjects, 2*3 items, W generated automatically,
#first parameter set to 0, no standard errors computed.

res <- LRSM(lrsmdat, mpoints = 2, groupvec = 1, sum0 = FALSE, se = FALSE)
res
}
\keyword{models}
