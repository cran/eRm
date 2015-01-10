\encoding{UTF-8}
\name{predict.ppar}
\alias{predict.ppar}

\title{Predict methods}
\description{Returns data matrix based on model probabilites. So far implemented for dichotomous models only.}
\usage{
\method{predict}{ppar}(object, cutpoint = "randomized", ...)
}

\arguments{
  \item{object}{Object of class \code{ppar} (from \code{person.parameter()}).}
  \item{cutpoint}{Either single integer value between 0 and 1 or \code{"randomized"} for randomized 0-1 assignment (see details)}
  \item{...}{Additional arguments ignored}
}

\details{
A randomized assignment implies that for each cell an additional random number is drawn. If the model probability is larger than this value, the person gets 1 on this particular item, if smaller, 0 is assigned. Alternatively, a numeric probability cutpoint can be assigned and the 0-1 scoring is carried out according to the same rule.
}
\value{
Returns data matrix based on model probabilities
}

\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{
    \code{\link{gofIRT.ppar}}
}

\examples{
#Model-based data matrix for RSM
res <- RM(raschdat2)
pres <- person.parameter(res)
predict(pres)
}
\keyword{models}
