\encoding{UTF-8}
\name{phi.range}
\alias{phi.range}
\title{Example User Function}
\description{Calculates the \eqn{R_\phi} statistic, i.e., the range of the inter-column correlations (\eqn{\phi}-coefficients) for a binary matrix.}
\usage{phi.range(mat)}
\arguments{\item{mat}{a binary matrix}}
\value{The range of the inter-column correlations}
\examples{
ctr <- rsctrl(burn_in = 10, n_eff = 5, step=10, seed = 123, tfixed = FALSE)
mat <- matrix(sample(c(0,1), 50, replace = TRUE), nr = 10)
rso <- rsampler(mat, ctr)
rso_st <- rstats(rso,phi.range)
print(unlist(rso_st))
}
\keyword{misc}
