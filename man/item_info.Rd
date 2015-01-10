\encoding{UTF-8}
\name{item_info}
\alias{item_info}
\alias{i_info}

\title{Calculate Item Information for 'eRm' objects
}
\description{
Calculates Samejima's (1969) information for all items
}
\usage{
item_info(ermobject, theta = seq(-5, 5, 0.01))

i_info(hvec, itembeta, theta)
}
\arguments{
  \item{ermobject}{An object of class 'eRm'.
  } 
  \item{theta}{Supporting or sampling points on the latent trait.
  }
  \item{hvec}{Number of categories of a single item.
  }
  \item{itembeta}{Cumulative item category parameters for a single item.
  }
}
\details{The function \code{item_info} calculates information of the
  whole set of items in the 'eRm' object. The function \code{i_info}
  does the same for a single item (and is called by \code{item_info}).  
}
\value{
  Returns a list (\code{i_info}) or a list of lists (where each element
  corresponds to an item, \code{item_info}) and contains
   \item{c.info}{Matrix of category information in columns for the
    different theta values in rows.}
   \item{i.info}{Vector of item information for the
    different theta values.}
}
\references{
Samejima, F. (1969) Estimation of latent ability using a response
pattern of graded scores. \emph{Psychometric Monographs}, \bold{17}.  
}
\author{Thomas Rusch} 
\seealso{
The function to calculate the test information, \code{\link{test_info}}
and the plot function \code{\link{plotINFO}}.
}
\examples{
res <- PCM(pcmdat)
info <- item_info(res)
plotINFO(res,type="item")
}
