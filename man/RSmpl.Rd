\encoding{UTF-8}
\name{RSmpl}
\alias{RSmpl}
\alias{RSmplext}
\title{Sample Objects}
\description{
  The objects of class \code{RSmpl} and \code{RSmplext} contain
  the original input matrix, the generated (encoded) random matrices, and
  some information about the sampling process.
}
\value{
   A list of class \code{RSmpl} or \code{RSmplext} with components
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
  \item{n_tot}{number of matrices in \code{outvec}.}
  \item{outvec}{vector of encoded random matrices}
  \item{ier}{error code (see below)}
}
\note{By default, all generated matrices plus
      the original matrix (in position 1) are contained in
      \code{outvec}, thus \code{n_tot = n_eff + 1}. If
      the original matrix is not in \code{outvec} then
      \code{n_tot = n_eff}.\cr\cr
      If \code{ier} is 0, no error was detected. Otherwise use
      the error function \code{rserror(ier)} to obtain some informations.\cr\cr
      For saving and loading objects
      of class \code{RSmpl} or \code{RSmplext}
      see the example in \code{\link{rsextrobj}}.
}
\section{Generation}{
      These classes of objects are returned from
      \code{rsampler} and \code{rsextrobj}.
}
\section{Methods}{
      Both classes have methods for the generic \code{summary}
      function.
}
\seealso{\code{\link{rsampler}}, \code{\link{rsextrobj}} }
\keyword{misc}
