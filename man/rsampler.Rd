\encoding{UTF-8}
\name{rsampler}
\alias{rsampler}
\title{Sampling Binary Matrices}
\description{%
  The function implements an \acronym{MCMC} algorithm for sampling of binary matrices with fixed margins complying to the Rasch model.
  Its stationary distribution is uniform.
  The algorithm also allows for square matrices with fixed diagonal.
}
\usage{
rsampler(inpmat, controls = rsctrl())
}
\arguments{
  \item{inpmat}{A binary (data) matrix with \eqn{n}{n} rows and \eqn{k}{k} columns.}
  \item{controls}{An object of class \code{\link{RSctr}}. If not specified, the default parameters as returned by function \code{\link{rsctrl}} are used.}
}
\details{
  \code{rsampler} is a wrapper function for a Fortran routine to generate binary random matrices based
   on an input matrix.
   On output the generated binary matrices are integer encoded. For further
   processing of the generated matrices use the function \code{\link{rstats}}.
}
\value{
   A list of class \code{\link{RSmpl}} with components
  \item{n}{number of rows of the input matrix}
  \item{k}{number of columns of the input matrix}
  \item{inpmat}{the input matrix}
  \item{tfixed}{\code{TRUE}, if diagonals of \code{inpmat} are fixed}
  \item{burn_in}{length of the burn in process}
  \item{n_eff}{number of generated matrices (effective matrices)}
  \item{step}{controls the number number of void matrices generated in the the burn in
                process and when effective matrices are generated (see note
                in \code{\link{rsctrl}}). }
  \item{seed}{starting value for the random number generator}
  \item{n_tot}{number of matrices in \code{outvec}, \code{n_tot = n_eff + 1}}
  \item{outvec}{vector of encoded random matrices}
  \item{ier}{error code}
}
\references{%
Verhelst, N. D. (2008). An Efficient \acronym{MCMC} Algorithm to Sample Binary Matrices with Fixed Marginals. \emph{Psychometrika, 73} (4)%
}
\author{Reinhold Hatzinger, Norman Verhelst}
\note{%
  An element of \code{outvec} is a four byte (or 32 bits) integer.
  The matrices to be output are stored bitwise (some bits are unused, since a integer is used for every row of a matrix).
  So the number of integers per row needed equals \eqn{(k+31)/32}{(k+31)/32} (integer division), which is one to four in the present implementation since the number of columns and rows must not exceed 128 and 4096, respectively.

  The summary method (\code{\link{summary.RSmpl}}) prints information on the content of the output object.
}
\seealso{\code{\link{rsctrl}}, \code{\link{rstats}} }
\examples{
data(xmpl)
ctr<-rsctrl(burn_in=10, n_eff=5, step=10, seed=0, tfixed=FALSE)
res<-rsampler(xmpl,ctr)
summary(res)
}
\keyword{misc}
