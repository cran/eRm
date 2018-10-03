\encoding{UTF-8}
\name{MLoef}
\alias{MLoef}
\alias{print.MLoef}
\alias{summary.MLoef}
\title{Martin-Löf's Likelihood-Ratio-Test}
\description{This Likelihood-Ratio-Test is based on item subgroup splitting.}
\usage{MLoef(robj, splitcr = "median")}
\arguments{
  \item{robj}{
    An object of class \code{'Rm'}.
  }
  \item{splitcr}{
    Split criterion to define the item groups.
    \code{"median"} and \code{"mean"} split items in two groups based on their items' raw scores.\cr%
    \code{splitcr} can also be a vector of length \eqn{k}{k} (where \eqn{k}{k} denotes the number of items) that takes two or more distinct values to define groups used for the Martin-Löf Test.
  }
}
\details{
  This function implements a generalization of the Martin-Löf test for polytomous items as proposed by Christensen, Bjørner, Kreiner & Petersen (2002), but does currently not allow for missing values.

  If the split criterion is \code{"median"} or \code{"mean"} and one or more items' raw scores are equal the median resp. mean, \code{MLoef} will assign those items to the lower raw score group.
  \code{summary.MLoef} gives detailed information about the allocation of all items.

  \code{summary} and \code{print} methods are available for objects of class \code{'MLoef'}.

  An \sQuote{exact} version of the Martin-Löf test for binary items is implemented in the \code{\link[eRm:NPtest]{NPtest}} function.
}
\value{
  \code{MLoef} returns an object of class \code{MLoef} containing:
    \item{LR}{LR-value}
    \item{df}{degrees of freedom}
    \item{p.value}{\emph{p}-value of the test}
    \item{fullModel}{the overall Rasch model}
    \item{subModels}{a list containing the submodels}
    \item{Lf}{log-likelihood of the full model}
    \item{Ls}{list of the sub models' log-likelihoods}
    \item{i.groups}{a list of the item groups}
    \item{splitcr}{submitted split criterion}
    \item{split.vector}{binary allocation of items to groups}
    \item{warning}{items equalling median or mean for the respective split criteria}
    \item{call}{the matched call}
}
\references{
Christensen, K. B., Bjørner, J. B., Kreiner S. & Petersen J. H. (2002). Testing unidimensionality in polytomous Rasch models. \emph{Psychometrika, (67)}4, 563--574.

Fischer, G. H., and Molenaar, I. (1995). \emph{Rasch Models -- Foundations, Recent Developements, and Applications.} Springer.

Rost, J. (2004). \emph{Lehrbuch Testtheorie -- Testkonstruktion.} Bern: Huber.
}
\author{Marco J. Maier, Reinhold Hatzinger}
\seealso{\code{\link{LRtest}}, \code{\link{Waldtest}}}
\examples{
# Martin-Löf-test on dichotomous Rasch model using "median" and a user-defined
# split vector. Note that group indicators can be of character and/or numeric.
splitvec <- c(1, 1, 1, "x", "x", "x", 0, 0, 1, 0)

res <- RM(raschdat1[,1:10])

MLoef.1 <- MLoef(res, splitcr = "median")
MLoef.2 <- MLoef(res, splitcr = splitvec)

MLoef.1

summary(MLoef.2)
}
\keyword{models}
