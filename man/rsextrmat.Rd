\encoding{UTF-8}
\name{rsextrmat}
\alias{rsextrmat}
\title{Extracting a Matrix}
\description{
  Convenience function to extract a matrix.
}
\usage{
rsextrmat(RSobj, mat.no = 1)
}
\arguments{
  \item{RSobj}{object as obtained from using \code{rsampler} or \code{rsextrobj}}
  \item{mat.no}{number of the matrix to extract from the sample object.}
}
\value{
   One of the matrices (either the original or a sampled matrix)
}
\seealso{\code{\link{rsampler}}, \code{\link{rsextrobj}},\code{\link{rstats}},}
\examples{
ctr <- rsctrl(burn_in = 10, n_eff = 3, step=10, seed = 0, tfixed = FALSE)
mat <- matrix(sample(c(0,1), 50, replace = TRUE), nr = 10)
all_m <- rsampler(mat, ctr)
summary(all_m)

# extract the third sampled matrix (here the fourth)
third_m <- rsextrmat(all_m, 4)
head(third_m)
}
\keyword{misc}
