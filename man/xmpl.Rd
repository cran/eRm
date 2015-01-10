\encoding{UTF-8}
\name{xmpl}
\alias{xmpl}
\alias{xmplbig}
\docType{data}
\title{Example Data}
\description{
  Ficitious data sets - matrices with binary responses
}
\usage{data(xmpl)}
\format{
  The format of \code{xmpl} is:\cr
  300 rows (referring to subjects) \cr
   30 columns (referring to items) \cr

  The format of \code{xmplbig} is:\cr
  4096 rows (referring to subjects) \cr
   128 columns (referring to items) \cr
  \code{xmplbig} has the maximum dimensions that the RaschSampler package
  can handle currently.
}
\examples{
data(xmpl)
print(head(xmpl))
}
\keyword{datasets}
