\encoding{UTF-8}
\name{sim.2pl}
%
\alias{sim.2pl}
%
%
%
\title{Simulation of 2-PL Data}
\description{This utility function returns a 0-1 matrix violating the parallel ICC assumption in the Rasch model.}
\usage{sim.2pl(persons, items, discrim = 0.25, seed = NULL, cutpoint = "randomized")}
%
\arguments{
  \item{persons}{Either a vector of person parameters or an integer indicating the number of persons (see details).}
  \item{items}{Either a vector of item parameters or an integer indicating the number of items (see details).}
  \item{discrim}{Standard deviation on the log scale.}
  \item{seed}{A seed for the random number generated can be set.}
  \item{cutpoint}{Either \code{"randomized"} for a randomized transformation of the model probability matrix into the model 0-1 matrix or an integer value between 0 and 1 (see details).}
}

\details{%
If \code{persons} and/or \code{items} (using single integers) are specified to determine the number of subjects or items, the corresponding parameter vector is drawn from N(0,1).
The \code{cutpoint} argument refers to the transformation of the theoretical probabilities into a 0-1 data matrix.
A randomized assingment implies that for each cell an additional random number is drawn.
If the model probability is larger than this value, the person gets 1 on this particular item, if smaller, 0 is assigned.
Alternatively, a numeric probability cutpoint can be assigned and the 0-1 scoring is carried out according to the same rule.

The \code{discrim} argument can be specified either as a vector of length \code{items} defining the item discrimination parameters in the 2-PL (e.g., \code{c(1,1,0.5,1,1.5)}), or as a single value.
In that case, the discrimination parameters are drawn from a lognormal distribution with \code{meanlog = 0}, where the specified value in \code{discrim} refers to the standard deviation on the log-scale.
The larger the values, the stronger the degree of Rasch violation.
Reasonable values are up to 0.5.
If 0, the data are Rasch homogeneous.
}
%
%
%
\references{
Su\'arez-Falc\'on, J. C., & Glas, C. A. W. (2003). Evaluation of global testing procedures for
   item fit to the Rasch model. British Journal of Mathematical and Statistical Society,
   56, 127-143.
}
%
\seealso{\code{\link{sim.rasch}}, \code{\link{sim.locdep}}, \code{\link{sim.xdim}}}
%
\examples{
#simulating 2-PL data
#500 persons, 10 items, sdlog = 0.30, randomized cutpoint
X <- sim.2pl(500, 10, discrim = 0.30)

#item and discrimination parameters from uniform distribution,
#cutpoint fixed
dpar <- runif(50, 0, 2)
ipar <- runif(50, -1.5, 1.5)
X <- sim.2pl(500, ipar, dpar, cutpoint = 0.5)
}
%
\keyword{models}
