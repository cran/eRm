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
\details{The Wald test only works if both subgroups have the same parameters. For instance,
for small samples in RSM or PCM the user must find an appropriate split such that
in each subgroup the same item-category parameters occur.
}
\value{
Returns an object of class \code{wald} containing:   
  \item{etapar1}{Eta parameters for first subgroup}
  \item{se1}{Standard errors for first subgroup}
  \item{etapar2}{Eta parameters for second subgroup}
  \item{se2}{Standard errors for second subgroup}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.

Fischer, G. H., and Scheiblechner, H. (1970). Algorithmen und Programme für das
probabilistische Testmodell von Rasch [Algorithms and programs for Rasch's
probabilistic test model]. Psychologische Beiträge, 12, 23-51.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{\code{\link{Waldtest}}}
\examples{

#Wald test for dichotomous Rasch model with median subject split
data(raschdat1)
res <- RM(raschdat1)
Waldtest(res)

#Wald test with user-defined subject split
splitvec <- sample(1:2,100,replace=TRUE)
Waldtest(res, splitcr = splitvec)

}

\keyword{models}
