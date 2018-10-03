\encoding{UTF-8}
\name{summary.RSctr}
\alias{summary.RSctr}
\title{Summary Method for Control Objects}
\description{
  Prints the current definitions for the sampling function.
}
\usage{
\method{summary}{RSctr}(object, ...)
}
\arguments{
  \item{object}{ object of class \code{RSctr} as obtained from \code{\link{rsctrl}} }
  \item{\dots}{ potential further arguments (ignored) }
}
\seealso{ \code{\link{rsctrl}} }
\examples{
   ctr <- rsctrl(n_eff = 1, seed = 123123123)  # specify controls
   summary(ctr)
}
\keyword{misc}
