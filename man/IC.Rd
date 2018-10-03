\encoding{UTF-8}
\name{IC}
\alias{IC}
\alias{IC.ppar}
\title{Information criteria}
\description{Computation of information criteria such as AIC, BIC, and cAIC based on
unconditional (joint), marginal, and conditional log-likelihood}
\usage{
\method{IC}{ppar}(object)
}

\arguments{
  \item{object}{Object of class \code{ppar} (from \code{person.parameter()}.}
}

\details{
The joint log-likelihood is established by summation of the logarithms of the estimated
solving probabilities. The marginal log-likelihood can be computed directly from the
conditional log-likelihood (see vignette for details).
}

\value{
The function \code{IC} returns an object of class \code{ICr} containing:
  \item{ICtable}{Matrix containing log-likelihood values, number of parameters, AIC, BIC, and
   cAIC for the joint, marginal, and conditional log-likelihood.}
}

\seealso{
    \code{\link{LRtest.Rm}}
}
\examples{
#IC's for Rasch model
res <- RM(raschdat2)             #Rasch model
pres <- person.parameter(res)    #Person parameters
IC(pres)

#IC's for RSM
res <- RSM(rsmdat)
pres <- person.parameter(res)
IC(pres)
}
\keyword{models}
