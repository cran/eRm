\name{person.parameter}
\alias{person.parameter}
\alias{person.parameter.eRm}
\alias{summary.ppar}
\alias{print.ppar}
\alias{plot.ppar}
\alias{logLik.ppar}
\alias{print.logLik.ppar}
\alias{IC}
\alias{IC.ppar}
\alias{confint.ppar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of Person Parameters}
\description{Maximum likelihood estimation of the person parameters with spline
interpolation for non-observed and 0/full responses. Extraction of information criteria such
as AIC, BIC, and cAIC based on unconditional log-likelihood.}
\usage{
\method{person.parameter}{eRm}(object)
\method{summary}{ppar}(object, ...)
\method{print}{ppar}(x, ...)
\method{summary}{ppar}(object, ...)
\method{plot}{ppar}(x, xlab = "Person Raw Scores", ylab = "Person Parameters (Theta)", main = NULL, ...)
\method{logLik}{ppar}(object, ...)
\method{IC}{ppar}(object)
\method{confint}{ppar}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{object}{Object of class \code{eRm} in \code{person.parameter} and object of class \code{ppar} in \code{IC}.}

Arguments for \code{print} and \code{plot} methods:
  \item{x}{Object of class \code{ppar}.}
  \item{xlab}{Label of the x-axis.}
  \item{ylab}{Label of the y-axis.}
  \item{main}{Title of the plot.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}

Arguments for \code{confint}:
  \item{parm}{Parameter specification (ignored).}
  \item{level}{Alpha-level.}
}
\details{If the data set contains missing values, person parameters are estimated
for each missing value subgroup.
}
\value{
The function \code{person.parameter} returns an object of class \code{ppar} containing:
  \item{loglik}{Log-likelihood of the collapsed data (for faster estimation persons with the same raw score are collapsed).}
  \item{npar}{Number of parameters.}
  \item{niter}{Number of iterations.}
  \item{thetapar}{Person parameter estimates.}
  \item{se.theta}{Standard errors of the person parameters.}
  \item{hessian}{Hessian matrix.}

The function \code{IC} returns an object of class \code{ICr} containing:
  \item{j.loglik}{Unconditional (joint) log-likelihood value.}
  \item{AIC}{AIC criterion.}
  \item{BIC}{BIC criterion.}
  \item{cAIC}{Consistent AIC criterion.}
}

\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations,
Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the eRm package in R. Psychology Science, 49, 26-43.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{
    \code{\link{itemfit.ppar}},\code{\link{personfit.ppar}}
}
\examples{

#Person parameter estimation of a rating scale model
data(rsmdat)
res <- RSM(rsmdat)
pres <- person.parameter(res)
print(pres)
summary(pres)
plot(pres)

#Person parameter estimation for a Rasch model with missing values
data(raschdat2)
res <- RM(raschdat2, se = FALSE)           #Rasch model without standard errors
pres <- person.parameter(res)
print(pres)                                #person parameters
summary(pres)
logLik(pres)                               #conditional log-lik value

#Extraction of IC's
IC(pres)

}
\keyword{models}
