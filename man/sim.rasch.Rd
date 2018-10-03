\encoding{UTF-8}
\name{sim.rasch}
\alias{sim.rasch}

\title{Simulation of Rasch homogeneous data}
\description{This utility function returns a 0-1 matrix which fits the Rasch model.
}
\usage{
sim.rasch(persons, items, seed = NULL, cutpoint = "randomized")
}

\arguments{
  \item{persons}{Either a vector of person parameters or an integer indicating the number of persons (see details)}
  \item{items}{Either a vector of item parameters or an integer indicating the number of items (see details)}
  \item{seed}{A seed for the random number generated can be set.}
  \item{cutpoint}{Either \code{"randomized"} for a randomized tranformation of the model probability matrix into the model 0-1 matrix or an integer value between 0 and 1 (see details)}
}

\details{If \code{persons} or \code{items} is an integer value, the corresponding parameter vector is drawn from N(0,1). The \code{cutpoint} argument refers to the transformation of the theoretical probabilities into a 0-1 data matrix. A randomized assingment implies that for each cell an additional random number is drawn. If the model probability is larger than this value, the person gets 1 on this particular item, if smaller, 0 is assigned. Alternatively, a numeric probability cutpoint can be assigned and the 0-1 scoring is carried out according to the same rule.
}

\references{
Su\'arez-Falc\'on, J. C., & Glas, C. A. W. (2003). Evaluation of global testing procedures for
   item fit to the Rasch model. British Journal of Mathematical and Statistical Society,
   56, 127-143.
}

\seealso{\code{\link{sim.xdim}}, \code{\link{sim.locdep}}, \code{\link{sim.2pl}}}
\examples{

#simulating Rasch homogenous data
#100 persons, 10 items, parameter drawn from N(0,1)
X <- sim.rasch(100, 10)

#person parameters drawn from uniform distribution, fixed cutpoint
ppar <- runif(100,-2,2)
X <- sim.rasch(ppar, 10, cutpoint = 0.5)
}

\keyword{models}

