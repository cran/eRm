\name{MLoef}
\alias{MLoef}
\alias{print.MLoef}
\alias{summary.MLoef}
\title{Computation of Martin-Loef's LR-Test}
\description{This LR-Test is based on item subgroup splitting.}
\usage{
MLoef(robj, splitcr = "median")
}
\arguments{
  \item{robj}{Object of class \code{Rm}.}
  \item{splitcr}{Split criterion to define the item groups.
    \code{"median"} and \code{"mean"} split items in two groups based on their
    items' raw scores. \code{splitcr} can also be a vector of length k (where k
    denotes the number of items) that takes two or more distinct values to
    define groups used for the Martin-Loef Test.}
}
\details{
  This function implements a generalization of the Martin-Loef test for polytomous
  items as proposed by Christensen, Bjorner, Kreiner & Petersen (2002), but does
  currently not allow for missing values.
    % The function can handle missing values, as long as every subject has at
    % least 2 valid responses in each group of items.

  If the split criterion is \code{"median"} or \code{"mean"} and one or more items'
  raw scores are equal the median resp. mean, \code{MLoef} will assign those items
  to the lower raw score group. \code{summary.MLoef} gives detailed information
  about the allocation of all items.

  \code{summary} and \code{print} methods are available for objects of class
  \code{MLoef}.

  An 'exact' version of the Martin-Loef test for binary items is implemented
  in the function \code{\link{NPtest}}.
}
\value{
  \code{MLoef} returns an object of class \code{MLoef} containing:
    \item{LR}{LR-value}
    \item{df}{degrees of freedom of the test statistic}
    \item{p.value}{p-value of the test}
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
Christensen, K. B., Bjorner, J. B., Kreiner S. & Petersen J. H. (2002). Testing unidimensionality in polytomous Rasch models. \emph{Psychometrika, (67)}4, 563--574.

Fischer, G. H., and Molenaar, I. (1995). \emph{Rasch Models -- Foundations, Recent Developements, and Applications.} Springer.

Rost, J. (2004). \emph{Lehrbuch Testtheorie -- Testkonstruktion.} Bern: Huber.
}
\author{Marco Maier, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{LRtest}}, \code{\link{Waldtest}}}
\examples{
# Martin-Loef-test on dichotomous Rasch model using "median" and a user-defined
# split vector. Note that group indicators can be of character and/or numeric.
splitvec <- c(1, 1, 1, "x", "x", "x", 0, 0, 1, 0)

res <- RM(raschdat1[,1:10])

MLoef.1 <- MLoef(res, splitcr = "median")
MLoef.2 <- MLoef(res, splitcr = splitvec)

MLoef.1

summary(MLoef.2)
}
\keyword{models}
