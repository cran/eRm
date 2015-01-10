\encoding{UTF-8}
\name{Analysis of Deviances}
\alias{anova.eRm}
\alias{print.eRm_anova}
\title{Analysis of Deviances for Rasch Models}
\description{Performs likelihood ratio tests against the model with the largest number of parameters.}
\usage{\method{anova}{eRm}(object, \dots)

\method{print}{eRm_anova}(x, \dots)}
\arguments{
  \item{object}{
    Gives the first object to be tested against others which follow, separated by commata.
  }
  \item{x}{
    An object of class \code{"eRm_anova"}.
  }
  \item{\dots}{
    Further models to test with \code{anova.eRm()}.
  }
}
\details{
The \code{anova} method is quite flexible and, as long the used data are identical, every model except the \acronym{LLRA} can be tested against each other.
Regardless of the order that models are specified, they will always be sorted by the number of parameters in decreasing order.
If \eqn{\geq3}{>= 3} models are passed to the method, all models will be tested against the first model (i.e., the one with the largest amount of parameters).
}
\value{
  \code{anova.eRm} returns a list object of class \code{eRm_anova} containing:
    \item{calls}{function calls of the different models (character).}
    \item{statistics}{the analysis of deviances table (columns are \code{LLs}: conditional log-likelihoods, \code{dev}: deviances, \code{npar}: number of parameters, \code{LR}: likelihood ratio statistics, \code{df}: degrees of freedom, \code{p}: \eqn{p}{p}-values).}
}
\section{Warning}{Although, there is a check for identical data matrices used, the models have to be nested for the likelihood ratio test to work.
You have to ensure that this is the case, otherwise results will be invalid.

\acronym{LLRA}s cannot be tested with other models (\acronym{RM}, \acronym{LLTM}, \acronym{RSM}, \ldots); for more information see \code{\link{anova.llra}}.}
%\references{xxx}
\author{Marco J. Maier}
\seealso{\code{\link{anova.llra}}, \code{\link{anova}}}
\examples{
### dichotomous data
dmod1 <- RM(lltmdat1)
dmod2 <- LLTM(lltmdat1, mpoints = 2)
anova(dmod1, dmod2)

### polytomous data
pmod1 <- RSM(rsmdat)
pmod2 <- PCM(rsmdat)
anova(pmod1, pmod2)

W <- cbind(rep(c(1,0), each=9), rep(c(0,1), each=9))
W
pmod3 <- LPCM(rsmdat, W)
anova(pmod3, pmod1, pmod2) # note that models are sorted by npar
}
