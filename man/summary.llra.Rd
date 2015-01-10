\encoding{UTF-8}
\name{summary.llra}
\alias{summary.llra}
\alias{print.summary.llra}

\title{Summarizing Linear Logistic Models with Relaxed Assumptions (LLRA)
}
\description{
\code{summary} method for class \code{"llra"}
}
\usage{
\method{summary}{llra}(object, level, ...)

\method{print}{summary.llra}(x, ...)
}
\arguments{
  \item{object}{an object of class "llra", typically result of a call to
    \code{\link{LLRA}}.
    }
  \item{x}{an object of class "summary.llra", usually, a result of a call
    to \code{summary.llra}.
  }
   \item{level}{The level of confidence for the confidence
   intervals. Default is 0.95.}
    \item{\dots}{further arguments passed to or from other methods.
    }
}
\details{
Objects of class \code{"summary.llra"} contain all parameters of interest plus the confidence intervals.

\code{print.summary.llra} rounds the values to 3 digits and displays
them nicely.
}
\value{
  The function \code{summary.lllra} computes and returns a list of
  summary statistics of the fitted LLRA given in object, reusing the
  components (list elements) \code{call}, \code{etapar},
  \code{iter}, \code{loglik}, \code{model}, \code{npar} and \code{se.etapar} from its argument, plus

  \item{ci}{The upper and lower confidence interval borders.}
}
\author{Thomas Rusch}

\seealso{
The model fitting function \code{\link{LLRA}}.
}
\examples{
##Example 6 from Hatzinger & Rusch (2009)
groups <- c(rep("TG",30),rep("CG",30))
llra1 <- LLRA(llradat3,mpoints=2,groups=groups)
summary(llra1)

\dontrun{
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.
ex2 <- LLRA(llraDat2[1:20],mpoints=4,llraDat2[21])
sumEx2 <- summary(ex2, level=0.95)

#print a summary
sumEx2

#get confidence intervals
sumEx2$ci}
}
