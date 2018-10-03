\encoding{UTF-8}
\name{sim.xdim}
\alias{sim.xdim}

\title{Simulation of multidimensional binary data}
\description{This utility function simulates a 0-1 matrix violating the
  unidimensionality assumption in the Rasch model.
}
\usage{
sim.xdim(persons, items, Sigma, weightmat, seed = NULL,
   cutpoint = "randomized")
}

\arguments{
  \item{persons}{Either a matrix (each column corresponds to a dimension) of person parameters or an integer indicating the number of persons (see details).}
  \item{items}{Either a vector of item parameters or an integer indicating the number of items (see details).}
  \item{Sigma}{A positive-definite symmetric matrix specifying the covariance matrix of the variables.}
  \item{weightmat}{Matrix for item-weights for each dimension (columns).}
  \item{seed}{A seed for the random number generated can be set.}
  \item{cutpoint}{Either \code{"randomized"} for a randomized tranformation of the model probability matrix into the model 0-1 matrix or an integer value between 0 and 1 (see details).}
}

\details{If \code{persons} is specified as matrix, \code{Sigma} is ignored. If \code{items} is
an integer value, the corresponding parameter vector is drawn from N(0,1).
The \code{cutpoint} argument refers to the transformation of the theoretical probabilities
into a 0-1 data matrix. A randomized assingment implies that for each cell an additional random
number is drawn. If the model probability is larger than this value, the person gets 1 on
this particular item, if smaller, 0 is assigned. Alternatively, a numeric probability
cutpoint can be assigned and the 0-1 scoring is carried out according to the same rule.

If \code{weightmat} is not specified, a random indicator matrix is generated where each item is a measurement
of only one dimension. For instance, the first row for a 3D-model could be (0,1,0) which means
that the first item measures the second dimension only. This corresponds to the between-item
multidimensional model presented by Adams et al. (1997).

\code{Sigma} reflects the VC-structure for the person parameters drawn from a multivariate
standard normal distribution. Thus, the diagonal elements are typically 1 and the lower the
covariances in the off-diagonal, the stronger the model violation.
}

\references{
Adams, R. J., Wilson, M., & Wang, W. C. (1997). The multidimensional random coefficients
   multinomial logit model. Applied Psychological Measurement, 21, 1-23.

Glas, C. A. W. (1992). A Rasch model with a multivariate distribution of ability. In M.
   Wilson (Ed.), Objective Measurement: Foundations, Recent Developments, and
   Applications (pp. 236-258). Norwood, NJ: Ablex.
}

\seealso{\code{\link{sim.rasch}}, \code{\link{sim.locdep}}, \code{\link{sim.2pl}}}
\examples{

# 500 persons, 10 items, 3 dimensions, random weights.
Sigma <- matrix(c(1, 0.01, 0.01, 0.01, 1, 0.01, 0.01, 0.01, 1), 3)
X <- sim.xdim(500, 10, Sigma)

#500 persons, 10 items, 2 dimensions, weights fixed to 0.5
itemvec <- runif(10, -2, 2)
Sigma <- matrix(c(1, 0.05, 0.05, 1), 2)
weights <- matrix(0.5, ncol = 2, nrow = 10)
X <- sim.xdim(500, itemvec, Sigma, weightmat = weights)

}
\keyword{models}

