\encoding{UTF-8}
\name{summary.RSmpl}
\alias{summary.RSmpl}
\alias{summary.RSmplext}
\title{Summary Methods for Sample Objects}

\description{
  Prints a summary list for sample objects of class \code{\link{RSmpl}}
  and \code{\link{RSmplext}}.
}
\usage{
\method{summary}{RSmpl}(object, ...)
\method{summary}{RSmplext}(object, ...)
}
\arguments{
  \item{object}{object as obtained from \code{rsampler} or \code{rsextrobj} }
  \item{\dots}{ potential further arguments (ignored) }
}
\details{
  Describes the status of an sample object.
}
\seealso{\code{\link{rsampler}}, \code{\link{rsextrobj}} }
\examples{
ctr <- rsctrl(burn_in = 10, n_eff = 3, step=10, seed = 0, tfixed = FALSE)
mat <- matrix(sample(c(0,1), 50, replace = TRUE), nr = 10)
all_m <- rsampler(mat, ctr)
summary(all_m)

some_m <- rsextrobj(all_m, 1, 2)
summary(some_m)
}
\keyword{misc}
