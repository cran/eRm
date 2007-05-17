\name{plotICC}
\alias{plotICC}
\alias{plotICC.Rm}
\alias{plotjointICC}
\alias{plotjointICC.dRm}

\title{ICC and Goodness-of-Fit Plots}
\description{Various plot functions for visualizing the item characteristic curves
}
\usage{
\method{plotICC}{Rm}(object, item.subset = "all", empirical = FALSE, xlim = c(-4, 4), ylim = c(0, 1), xlab = "Latent Dimension", ylab = "Probability to Solve", ...)
\method{plotjointICC}{dRm}(object, item.subset = "all", legend = TRUE, xlim = c(-4, 4), ylim = c(0, 1), xlab = "Latent Dimension", ylab = "Probability to Solve", lty = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{Rm} or \code{dRm}}
  \item{item.subset}{Subset of items to be plotted. Either a numeric vector indicating
  the column in \code{X} or a character vector indiciating the column name. 
  If \code{"all"}, all items are plotted.}
  \item{empirical}{Option for plotting the empirical ICCs for objects of class \code{dRm}}
  \item{legend}{If \code{TRUE}, legend is provided, otherwise the ICCs are labeled.}
  \item{xlab}{Label of the x-axis.}
  \item{ylab}{Label of the y-axis.}
  \item{xlim}{Range of person parameters.}
  \item{ylim}{Range for probability to solve.}
  \item{lty}{Line type.}
  \item{...}{Additional plot parameters.}
}
\details{}
\value{}
\references{
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{}
\seealso{\code{\link{plotGOF}}}
\examples{

# Rating scale model, ICC plot for all items
data(rsmdat)
res <- RSM(rsmdat)
plotICC(res)

# Rasch model with empirical ICCs
newdata <- matrix(sample(0:1, 1000, replace = TRUE), ncol = 5)
res <- RM(newdata)
plotICC(res, empirical = TRUE)

# Joint ICC plot for items 2, 6, 8, and 15 for a Rasch model
data(raschdat1)
res <- RM(raschdat1)
plotjointICC(res, item.subset = c(2,6,8,15))

}
\keyword{models}
