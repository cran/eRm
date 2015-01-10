\encoding{UTF-8}
\name{plotGR}
\alias{plotGR}

\title{Plot Treatment or Covariate Effects for LLRA
}
\description{
Plots treatment or covariate group effects over time.  
}
\usage{
plotGR(object, ...)
}
\arguments{
  \item{object}{an object of class "llra".
  }
  \item{\dots}{Additional parameters to be passed to and from other
    methods.
  }
}
\details{
The plot is a lattice plot with each panel corresponding to an item. The
effects are plotted for each groups (including baseline) over the
different time points. The groups are given the same names as for the
parameter estimates (derived from groupvec).

Please note that all effects are to be interpreted relative to the
baseline.

Currently, this function only works for a full item x treatment x
timepoints LLRA. Collapsed effects will not be displayed properly. 
}
\author{Thomas Rusch}
\seealso{
The plot method for trend effects \code{\link{plotTR}}.
}
\section{Warning:}{
  Objects of class \code{"llra"} that contain estimates from a collapsed
  data matrix will not be displayed correctly.
}
\examples{
##Example 6 from Hatzinger & Rusch (2009)
groups <- c(rep("TG",30),rep("CG",30))
llra1 <- LLRA(llradat3,mpoints=2,groups=groups)
summary(llra1)
plotGR(llra1)

\dontrun{  
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.
ex2 <- LLRA(llraDat2[1:20],mpoints=4,groups=llraDat2[21])
plotGR(ex2)}
}
