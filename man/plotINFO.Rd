\name{plotINFO}
\alias{plotINFO}

\title{Plot Information For 'eRm' objects
}
\description{
Calculates and plots the individual or summed item information by Samejima (1969)
}
\usage{
plotINFO(ermobject, type="both", theta=seq(-5,5,0.01),...)
}
\arguments{
  \item{ermobject}{An object of class 'eRm'.
  }
  \item{type}{A string denoting the type of information to be
    plotted. Currently supports 'item', 'test' and 'both' (default).
  }
  \item{theta}{Supporting or sampling points on the latent trait.
  }
  \item{...}{Further arguments.}
}

\references{
Samejima, F. (1969) Estimation of latent ability using a response
pattern of graded scores. \emph{Psychometric Monographs}, \bold{17}.  
}
\author{
Thomas Rusch
} 
\seealso{
The function to calculate the item or test information,
\code{\link{item_info}} and \code{\link{test_info}}.
}
\examples{
res <- PCM(pcmdat)
plotINFO(res)
}



