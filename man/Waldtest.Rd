\name{Waldtest}
\alias{Waldtest}
\alias{Waldtest.Rm}
\alias{print.wald}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Item-Specific Wald Test}
\description{Performs a Wald test on item-level by splitting subjects into subgroups.
}
\usage{
\method{Waldtest}{Rm}(object, splitcr = "median")
\method{print}{wald}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{RM}.}
  \item{splitcr}{Split criterion for subject raw score splitting. \code{median} 
  uses the median as split criterion, \code{mean} performs a mean-split. 
  Optionally \code{splitcr} can also be a dichotomous vector which assigns each person to a 
  certain subgroup (e.g., following an external criterion). }
  \item{x}{Object of class \code{wald}.}
  \item{...}{Further arguments passed to or from other methods. They are ignored in this function.}
}
\details{Items are eliminated if they not have the same number of categories in each subgroup. 
To avoid this problem, for RSM and PCM it is considered to use a random or another user-defined split.
}
\value{
Returns an object of class \code{wald} containing:   
  \item{coef.table}{Data frame with test statistics, z- and p-values.}
  \item{betapar1}{Beta parameters for first subgroup}
  \item{se.beta1}{Standard errors for first subgroup}
  \item{betapar2}{Beta parameters for second subgroup}
  \item{se.beta2}{Standard errors for second subgroup}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.

Fischer, G. H., and Scheiblechner, H. (1970). Algorithmen und Programme fuer das
probabilistische Testmodell von Rasch [Algorithms and programs for Rasch's
probabilistic test model]. Psychologische Beitraege, 12, 23-51.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{\code{\link{LRtest}}}
\examples{

#Wald test for Rasch model with user-defined subject split
data(raschdat2)
res <- RM(raschdat2)
splitvec <- sample(1:2,25,replace=TRUE)
Waldtest(res, splitcr = splitvec)

#Wald test for RSM eliminates 4 items (with median split)
data(rsmdat)
res <- RSM(rsmdat)
Waldtest(res)

}

\keyword{models}
