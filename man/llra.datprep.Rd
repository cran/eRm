\encoding{UTF-8}
\name{llra.datprep}
\alias{llra.datprep}

\title{Prepare Data Set for LLRA Analysis
}
\description{
Converts wide data matrix in long format, sorts subjects according to
groups and builds assigment vector.
}
\usage{
llra.datprep(X, mpoints, groups, baseline)
}
\arguments{
  \item{X}{Data matrix as described in Hatzinger and Rusch (2009). It
    must be of wide format, e.g. for each person all item answers are
    written in columns for t1, t2, t3 etc. Hence each row corresponds to
    all observations for a single person.
    Missing values are not allowed.  
  }
    \item{mpoints}{The number of time points.
  }
  \item{groups}{Vector, matrix or data frame with subject/treatment
    covariates.
  }
  \item{baseline}{An optional vector with the baseline values for the
    columns in group.}  
}
\details{The function converts a data matrix from wide to long fromat as
  needed for LLRA. Additionally it sorts the subjects according to the
  different treatment/covariate groups. The group with the lowest
  (alpha-)numerical value will be the
  baseline.

  Treatment and covariate groups are either defined by a vector, or by a
  matrix or data frame. The latter will be combined to a vector of
  groups corresponding to a combination of each factor level per column
  with the factor levels of the other column. The (constructed or
  passed) vector will then be used to create the assignment vector.  
}
\value{
Returns a list with the components
\item{X}{Data matrix in long format with subjects sorted by groups.}
\item{assign.vec}{The assignment vector.}
\item{grp_n}{A vector of the number of subjects in each group.}
}
\author{Reinhold Hatzinger}
\seealso{
The function that uses this is \code{\link{LLRA}}. The values from
  \code{llra.datprep} can be passed to \code{\link{build_W}}.
}
\examples{
    # example 3 items, 3 timepoints, n=10, 2x2 treatments
    dat<-sim.rasch(10,9)
    tr1<-sample(c("a","b"),10,r=TRUE)
    tr2<-sample(c("x","y"),10,r=TRUE)

    # one treatment
    res<-llra.datprep(dat,mpoints=3,groups=tr1)
    res<-llra.datprep(dat,mpoints=3,groups=tr1,baseline="b") 

    # two treatments
    res<-llra.datprep(dat,mpoints=3,groups=cbind(tr1,tr2))
    res<-llra.datprep(dat,mpoints=3,groups=cbind(tr1,tr2),baseline=c("b","x")) 

    # two treatments - data frame
    tr.dfr<-data.frame(tr1, tr2)
    res<-llra.datprep(dat,mpoints=3,groups=tr.dfr) 
}
