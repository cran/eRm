\encoding{UTF-8}
\name{rsextrobj}
\alias{rsextrobj}
\title{Extracting Encoded Sample Matrices}
\description{
  Utility function to extract some of the generated matrices, still in encoded form.
}
\usage{
rsextrobj(RSobj, start = 1, end = 8192)
}
\arguments{
  \item{RSobj}{object as obtained from using \code{rsampler}}
  \item{start}{number of the matrix to start with. When specifying 1
              (the default value) the original input matrix is
               included in the output object.
              }
  \item{end}{last matrix to be extracted. If \code{end}
             is not specified, all matrices from \code{RSobj}
             are extracted (the maximal value is 8192, see
             \code{rsctrl}). If \code{end} is larger than
             the number of matrices stored in \code{RSobj},
             \code{end} is set to the highest possible value
             (i.e., \code{n_tot}).
            }
}
\value{
   A list of class \code{\link{RSmpl}} with components
  \item{n}{number of rows of the input matrix}
  \item{k}{number of columns of the input matrix}
  \item{inpmat}{the input matrix}
  \item{tfixed}{\code{TRUE}, if diagonals of \code{inpmat} are fixed}
  \item{burn_in}{length of the burn in process}
  \item{n_eff}{number of generated matrices (effective matrices)}
  \item{step}{controls the number number of void matrices generated in the burn in
              process and when effective matrices are generated (see note
              in \code{\link{rsctrl}}). }
  \item{seed}{starting value for the random number generator}
  \item{n_tot}{number of matrices in \code{outvec}.}
  \item{outvec}{vector of encoded random matrices}
  \item{ier}{error code}
}
\note{By default, all generated matrices plus
      the original matrix (in position 1) are contained in
      \code{outvec}, thus \code{n_tot = n_eff + 1}. If
      the original matrix is not in \code{outvec} then
      \code{n_tot = n_eff}.\cr
      For saving and loading objects
      of class \code{RSobj} see the example below.

      For extracting a decoded (directly usable) matrix use \code{\link{rsextrmat}}.
}
\seealso{\code{\link{rsampler}}, \code{\link{rsextrmat}} }
\examples{
ctr <- rsctrl(burn_in = 10, n_eff = 3, step=10, seed = 0, tfixed = FALSE)
mat <- matrix(sample(c(0,1), 50, replace = TRUE), nr = 10)
all_m <- rsampler(mat, ctr)
summary(all_m)

some_m <- rsextrobj(all_m, 1, 2)
summary(some_m)

\dontrun{
save(some_m, file = "some.RSobj.RData")
rm(some_m)
ls()

load("some.RSobj.RData")
summary(some_m)}
}
\keyword{misc}
