\encoding{UTF-8}
\name{plotINFO}
\alias{plotINFO}
\title{Plot Information For \code{'eRm'} objects
}
\description{Calculates and plots the individual item-category information
  (type='category'), item information (type='item') or test/scale information (i.e., summed item information, type='scale' or 'test') ) as defined by Samejima (1969)}
\usage{
plotINFO(ermobject, type = "both", theta = seq(-6, 6, length.out = 1001L), 
legpos = "topright", ...)
}
\arguments{
  \item{ermobject}{An object of class \code{'eRm'}.}
  \item{type}{A string denoting the type of information to be
    plotted. Currently supports \code{"category"}, \code{"item"},
    \code{"test"}, \code{"scale"} and \code{"both"} (which gives item
    and scale information; default).}
  \item{theta}{Supporting or sampling points on the latent trait.}
  \item{legpos}{Defines the positioning of the legend, as in \code{\link[eRm]{plotICC}}.}
  \item{...}{%
    Further arguments.
    \code{xlab} sets the label of the \eqn{x}{x} axis.
    \code{ylabI} and \code{ylabT} control the labeling of the item or test information plot.
    \code{mainI} and \code{mainT} set the titles for item/test information plots.
  }
}
\references{Samejima, F. (1969) Estimation of latent ability using a response pattern of graded scores. \emph{Psychometric Monographs}, \bold{17}.}
\author{Thomas Rusch}
\seealso{The function to calculate the item or test information, \code{\link[eRm]{item_info}} and \code{\link[eRm]{test_info}}.}
\examples{
res <- PCM(pcmdat)
plotINFO(res)
}
