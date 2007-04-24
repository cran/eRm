\name{LRtest}
\alias{LRtest.Rm}
\alias{LRtest}
\alias{print.LR}
\alias{summary.LR}
\alias{plotGOF}
\alias{plotGOF.LR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio model test}
\description{Computation of Andersen's LR-test. 
}
\usage{
\method{LRtest}{Rm}(object, splitcr = "median", alpha = 0.05, se = TRUE)
\method{plotGOF}{LR}(x, beta.subset = "all", xlab = "Beta Group 1", 
ylab = "Beta Group 2", ylim = c(-3, 3), xlim = c(-3, 3), type = "p", ...)
\method{print}{LR}(x,...)
\method{summary}{LR}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  
The parameters for \code{LRtest} which returns an object of class \code{LR} are: 
  \item{object}{Object of class \code{Rm}.}
  \item{splitcr}{Split criterion for subject raw score splitting. \code{all.r} corresponds to a 
  full raw score split, \code{median} uses the median as split criterion, \code{mean} performs a mean-split. 
  Optionally \code{splitcr} can also be a vector which assigns each person to a 
  certain subgroup (e.g., following an external criterion). }
  \item{alpha}{Significance level for Chi-square statistic.}
  \item{se}{If \code{TRUE}, standard errors are computed}
  
The parameters for \code{plotGOF} are:   
  \item{x}{Object of class \code{LR} for visualizing the fit of single items.}
  \item{beta.subset}{Either a numeric vector or a character vector indicating the columns of the 
  data matrix which should be plotted. If "all", all items are plotted.}
  \item{xlab}{Label of the x-axis.}
  \item{ylab}{Label of the y-axis.}
  \item{xlim}{Parameter range for group 1.}
  \item{ylim}{Parameter range for group 2.}
  \item{type}{Type of plot desired. By default points are plotted.}
  \item{...}{Additional plot parameters.}
  
}
\details{}
\value{
  \code{LRtest} returns an object of class \code{LR} containing:
  \item{LR}{LR-value.}
  \item{df}{Degrees of freedom of the test statistic.}
  \item{Chisq}{Chi-square value with corresponding df and alpha.}
  \item{pvalue}{P-value of the test.}
  \item{likgroup}{Log-likelihood values for the subgroups}
  \item{betalist}{List of beta parameters for the subgroups.}
}
\references{
Andersen, E. B. (1973). A goodness of fit test for the Rasch model. Psychometrika, 38,
123-140.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{\code{\link{Waldtest}}}
\examples{

# LR-test on dichotomous Rasch model with median split
data(raschdat1)
res <- RM(raschdat1)
lrres <- LRtest(res)
print(lrres)
summary(lrres)
plotGOF(lrres)

# LR-test on dichotomous Rasch model with user-specified split
splitvec <- sample(1:3, 100, replace = TRUE)          #3 random subgroups
lrres <- LRtest(res, splitcr = splitvec)
print(lrres)
summary(lrres)
plotGOF(lrres)

}

\keyword{models}
