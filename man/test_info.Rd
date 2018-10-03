\encoding{UTF-8}
\name{test_info}
\alias{test_info}

\title{Calculate Test Information For \code{eRm} objects
}
\description{
Calculates the information of a test or a scale as the sum of Samejima's (1969) information for all items.
}
\usage{
test_info(ermobject, theta=seq(-5,5,0.01))
}
\arguments{
  \item{ermobject}{An object of class \code{'eRm'}.
  } 
  \item{theta}{Supporting or sampling points on the latent trait.
  }
}
\details{The function \code{test_info} calculates the test or scale information of the
  whole set of items in the \code{'eRm'} object. 
}
\value{
  Returns the vector of test information for all values of theta.  
}
\references{
Samejima, F. (1969) Estimation of latent ability using a response
pattern of graded scores. \emph{Psychometric Monographs}, \bold{17}.  
}
\author{Thomas Rusch} 
\seealso{
The function to calculate the item information, \code{\link{item_info}}
and the plot function \code{\link{plotINFO}}.

}
\examples{
res <- PCM(pcmdat)
tinfo <- test_info(res)
plotINFO(res, type="test")
}



