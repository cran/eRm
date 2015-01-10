\encoding{UTF-8}
\name{rstats}
\alias{rstats}
\title{Calculating Statistics for the Sampled Matrices}
\description{
   This function is used to calculate user defined statistics for the
   (original and) sampled matrices. A user defined function has to
   be provided.
}
\usage{
rstats(RSobj, userfunc, ...)
}
\arguments{
  \item{RSobj}{object as obtained from using \code{\link{rsampler}}
               or \code{\link{rsextrobj}} }
  \item{userfunc}{a user defined function which performs operations
     on the (original and) sampled matrices. The first argument in the definition
     of the user function must be an object of type matrix.}
  \item{...}{further arguments, that are passed to the user function}
}
\value{
      A list of objects as specified in the user supplied function
}
\note{The encoded matrices that are contained in the
      input object \code{RSobj} are decoded and passed to the user function in turn.
      If \code{RSobj} is not an object obtained from either \code{\link{rsampler}}
      or \code{\link{rsextrobj}} or
      no user function is specified an error message is printed.
      A simple user function, \code{\link{phi.range}}, is included in
      the RaschSampler package for demonstration purposes.\cr

      \code{rstats} can be used to obtain the 0/1 values for any
      of the sampled matrices (see second example below). Please note,
      that the output from the user function is stored in a list where
      the number of components corresponds to the number of matrices passed
      to the user function (see third example).
}
\seealso{\code{\link{rsampler}}, \code{\link{rsextrobj}} }
\examples{
ctr <- rsctrl(burn_in = 10, n_eff = 5, step=10, seed = 12345678, tfixed = FALSE)
mat <- matrix(sample(c(0,1), 50, replace = TRUE), nr = 10)
rso <- rsampler(mat, ctr)
rso_st <- rstats(rso,phi.range)
unlist(rso_st)

# extract the third generated matrix
# (here, the first is the input matrix)
# and decode it into rsmat

rso2 <- rsextrobj(rso,4,4)
summary(rso2)
rsmat <- rstats(rso2, function(x) matrix(x, nr = rso2$n))
print(rsmat[[1]])

# extract only the first r rows of the third generated matrix

mat<-function(x, nr = nr, r = 3){
  m <- matrix(x, nr = nr)
  m[1:r,]
}
rsmat2 <- rstats(rso2, mat, nr=rso$n, r = 3)
print(rsmat2[[1]])

# apply a user function to the decoded object
print(phi.range(rsmat[[1]]))
}
\keyword{misc}
