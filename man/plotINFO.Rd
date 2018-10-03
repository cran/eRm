\encoding{UTF-8}
\name{plotINFO}
\alias{plotINFO}
\title{Plot Information For \code{'eRm'} objects
}
\description{Calculates and plots the individual or summed item information by Samejima (1969)}
\usage{plotINFO(ermobject, type = "both", theta = seq(-6, 6, length.out = 1001L), \dots)}
\arguments{
  \item{ermobject}{An object of class \code{'eRm'}.}
  \item{type}{A string denoting the type of information to be plotted. Currently supports \code{"item"}, \code{"test"} and \code{"both"} (default).}
  \item{theta}{Supporting or sampling points on the latent trait.}
  \item{...}{%
    Further arguments.
    \code{xlab} sets the label of the \eqn{x}{x} axis.
    \code{ylabI} and \code{ylabT} control the labeling of the item or test information plot.
    \code{mainI} and \code{mainT} set the titles for item/test information plots.
    \code{legpos} defines the positioning of the legend, as in \code{\link[eRm]{plotICC}}.
  }
}
\references{Samejima, F. (1969) Estimation of latent ability using a response pattern of graded scores. \emph{Psychometric Monographs}, \bold{17}.}
\author{Thomas Rusch}
\seealso{The function to calculate the item or test information, \code{\link[eRm]{item_info}} and \code{\link[eRm]{test_info}}.}
\examples{
res <- PCM(pcmdat)
plotINFO(res)
}
