\name{person.parameter}
\alias{person.parameter}
\alias{person.parameter.eRm}
\alias{summary.ppar}
\alias{print.ppar}
\alias{plot.ppar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of Person Parameters}
\description{Maximum likelihood estimation of the person parameters.
}
\usage{
\method{person.parameter}{eRm}(object, se = TRUE, splineInt = TRUE)
\method{summary}{ppar}(object, ...)
\method{print}{ppar}(x, ...)
\method{plot}{ppar}(x, xlab = "Person Raw Scores", ylab = "Person Parameters (Theta)", main = NULL, ...)
}

\arguments{
  \item{object}{Object of class \code{eRm}.}
  \item{se}{If \code{TRUE}, the standard errors are estimated.}
  \item{splineInt}{\code{TRUE} if the missing person raw scores (0 and full raw score included)
  are interpolated using a cubic spline regression.} 
  
Arguments for \code{print} and \code{plot} method:
  \item{x}{Object of class \code{ppar}.}
  \item{xlab}{Label of the x-axis.}
  \item{ylab}{Label of the y-axis.}
  \item{main}{Title of the plot.}
  \item{...}{Further arguments to be passed to or from other methods. They are ignored in this function.}
}
\details{If the data set contains missing values, person parameters are estimated 
for each missing value subgroup.}

}
\value{
The function \code{person.parameter} returns an object of class \code{ppar} containing:
  \item{loglik}{Log-likelihood.}
  \item{npar}{Number of parameters.}
  \item{niter}{Number of iterations.}
  \item{thetapar}{Person parameter estimates.}
  \item{se.theta}{Standard errors of the person parameters.}  
}

\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the eRm package in R. Psychology Science, 49, 26-43.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{\code{\link{itemfit.ppar}},\code{\link{personfit.ppar}}}
}
\examples{

#Person parameter estimation of a rating scale model
data(rsmdat)
res <- RSM(rsmdat)
pres <- person.parameter(res)
print(pres)
summary(pres)
plot(pres)

#Person parameter estimation for a Rasch model with missing values and
#witout spline interpolation.
data(raschdat1)
raschNA <- raschdat1
raschNA[1:40,1] <- NA                      #first item not presented
raschNA[41:100,2] <- NA                    #second item not presented
res <- RM(raschNA, se = FALSE)             #Rasch model without standard errors
pres <- person.parameter(res, splineInt = FALSE)
print(pres)                                #person parameteres round to 5 digits
summary(pres)

}
\keyword{models}
