\encoding{UTF-8}
\name{gofIRT}
\alias{gofIRT}
\alias{gofIRT.ppar}
\alias{summary.gof}
\alias{print.gof}



\title{Various model tests and fit indices}

\description{This function computes various model tests and fit indices for objects of class \code{ppar}: Collapsed deviance, Casewise deviance, Rost's LR-test, Hosmer-Lemeshow test, R-Squared measures, confusion matrix, ROC analysis.
}

\usage{
\method{gofIRT}{ppar}(object, groups.hl = 10, cutpoint = 0.5)
}

\arguments{
  \item{object}{Object of class \code{ppar} (from \code{person.parameter()}).}
  \item{groups.hl}{Number of groups for Hosmer-Lemeshow test (see details).} 
  \item{cutpoint}{Integer between 0 and 1 for computing the 0-1 model matrix from the estimated probabilities}


}


\details{So far this test statistics are implemented only for dichotomous models without NA's. The Hosmer-Lemeshow test is computed by splitting the response vector into percentiles, e.g. \code{groups.hl = 10} corresponds to decile splitting. 
}
\value{
The function \code{gofIRT} returns an object of class \code{gof} containing:

  \item{test.table}{Ouput for model tests.}
  \item{R2}{List with R-squared measures.}
  \item{classifier}{Confusion matrix, accuracy, sensitivity, specificity.}
  \item{AUC}{Area under ROC curve.}
  \item{Gini}{Gini coefficient.}
  \item{ROC}{FPR and TPR for different cutpoints.}
  \item{opt.cut}{Optimal cutpoint determined by ROC analysis.}
  \item{predobj}{Prediction output from ROC analysis (\code{ROCR} package)}
} 

\references{
Mair, P., Reise, S. P., and Bentler, P. M. (2008). IRT goodness-of-fit using approaches from logistic regression. UCLA Statistics Preprint Series. 
}

\seealso{
    \code{\link{itemfit.ppar}},\code{\link{personfit.ppar}},\code{\link{LRtest}}
}
\examples{
#Goodness-of-fit for a Rasch model
res <- RM(raschdat1)
pres <- person.parameter(res)
gof.res <- gofIRT(pres)
gof.res
summary(gof.res)
}
\keyword{models}
