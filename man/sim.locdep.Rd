\encoding{UTF-8}
\name{sim.locdep}
\alias{sim.locdep}

\title{Simulation locally dependent items}
\description{This utility function returns a 0-1 matrix violating the
local independence assumption.
}
\usage{
sim.locdep(persons, items, it.cor = 0.25, seed = NULL,
   cutpoint = "randomized")
}

\arguments{
  \item{persons}{Either a vector of person parameters or an integer indicating the number of persons (see details).}
  \item{items}{Either a vector of item parameters or an integer indicating the number of items (see details).}
  \item{it.cor}{Either a single correlation value between 0 and 1 or a positive semi-definite VC matrix.}
  \item{seed}{A seed for the random number generated can be set.}
  \item{cutpoint}{Either \code{"randomized"} for a randomized tranformation of the model probability matrix into the model 0-1 matrix or an integer value between 0 and 1 (see details).}
}

\details{If \code{persons} or \code{items} is an integer value, the corresponding parameter vector
is drawn from N(0,1). The \code{cutpoint} argument refers to the transformation of the theoretical
probabilities into a 0-1 data matrix. A randomized assingment implies that for each cell an
additional random number is drawn. If the model probability is larger than this value,
the person gets 1 on this particular item, if smaller, 0 is assigned. Alternatively, a numeric probability cutpoint can be assigned and the 0-1 scoring is carried out according to the same rule.

The argument \code{it.cor} reflects the pair-wise inter-item correlation. If this should be constant
across the items, a single value between 0 (i.e. Rasch model) and 1 (strong violation) can be specified.
Alternatively, a symmetric VC-matrix of dimension number of items can be defined.

}

\references{
Jannarone, R. J. (1986). Conjunctive item response theory kernels. Psychometrika, 51,
357-373.

Su\'arez-Falc\'on, J. C., & Glas, C. A. W. (2003). Evaluation of global testing procedures for
   item fit to the Rasch model. British Journal of Mathematical and Statistical Society,
   56, 127-143.
}

\seealso{\code{\link{sim.rasch}}, \code{\link{sim.2pl}}, \code{\link{sim.xdim}}}
\examples{

#simulating locally-dependent data
#500 persons, 10 items, inter-item correlation of 0.5
X <- sim.locdep(500, 10, it.cor = 0.5)

#500 persons, 4 items, correlation matrix specified
sigma <- matrix(c(1,0.2,0.2,0.3,0.2,1,0.4,0.1,0.2,0.4,1,0.8,0.3,0.1,0.8,1),
   ncol = 4)
X <- sim.locdep(500, 4, it.cor = sigma)
}

\keyword{models}
