\name{LRtest}
\alias{LRtest.Rm}
\alias{LRtest}
\alias{print.LR}
\alias{summary.LR}
\alias{plotGOF}
\alias{plotGOF.LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of Andersen's LR-test.}
\description{This LR-test is based on subject subgroup splitting.
}
\usage{
\method{LRtest}{Rm}(object, splitcr = "median", se = FALSE)
\method{plotGOF}{LR}(x, beta.subset = "all", xlab = "Beta Group 1",
ylab = "Beta Group 2", tlab = "item", ylim = c(-3, 3), xlim = c(-3, 3), type = "p", pos = "4", ...)
\method{print}{LR}(x,...)
\method{summary}{LR}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{Rm}.}
  \item{splitcr}{Split criterion for subject raw score splitting. \code{all.r} corresponds to a
  full raw score split, \code{median} uses the median as split criterion, \code{mean} performs a mean-split.
  Optionally \code{splitcr} can also be a vector which assigns each person to a
  certain subgroup (e.g., following an external criterion). }
  \item{se}{If \code{TRUE} standard errors for beta's are computed.}

Arguments for \code{plotGOF}:
  \item{x}{Object of class \code{LR} for visualizing the fit of single items.}
  \item{beta.subset}{If \code{"all"}, all items are plotted. Otherwise numeric subset vector can be specified.}
  \item{xlab}{Label on x-axis.}
  \item{ylab}{Label on y-axis.}
  \item{tlab}{Specification of item labels: \code{"item"} prints the item names, \code{"number"} gives integers
       corresponding to order of the beta parameters, if \code{"none"} no labels are printed.
       \code{"identify"} allows for an interactive labelling. Initially no labels are printed, after clicking
       close to an item point the corresponding label is added. The identification process is terminated by clicking
       the second button and selecting 'Stop' from the menu, or from the 'Stop' menu on the graphics window.
       For more information see \code{\link{identify}}.
       }
  \item{xlim}{Limits on x-axis.}
  \item{ylim}{Limits on y-axis.}
  \item{type}{Plotting type.(see \code{\link{plot}})}
  \item{pos}{Position of the item label (see \code{\link{text}})}
  \item{...}{Additional graphical parameters.}
}
\details{If the data set contains missing values and \code{mean} or \code{median} is specified as splitcriterion,
         means or medians are calculated for each missing value subgroup and consequently used for raw score splitting.
}
\value{
\code{LRtest} returns an object of class \code{LR} containing:
  \item{LR}{LR-value.}
  \item{df}{Degrees of freedom of the test statistic.}
  \item{Chisq}{Chi-square value with corresponding df.}
  \item{pvalue}{P-value of the test.}
  \item{likgroup}{Log-likelihood values for the subgroups}
  \item{betalist}{List of beta parameters for the subgroups.}
  \item{selist}{List of standard errors of beta's.}
  \item{etalist}{List of eta parameters for the subgroups.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the eRm package in R. Psychology Science, 49, 26-43.
}
\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{Waldtest}}}
\examples{

# LR-test on dichotomous Rasch model with user-defined split
splitvec <- sample(1:3, 100, replace = TRUE)
data(raschdat1)
res <- RM(raschdat1)
lrres <- LRtest(res, splitcr = splitvec)
lrres
summary(lrres)
plotGOF(lrres)
}

\keyword{models}
