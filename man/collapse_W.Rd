\encoding{UTF-8}
\name{collapse_W}
\alias{collapse_W}

\title{
Convenient Collapsing of LLRA Design Matrix 
}
\description{
Collapses columns of a design matrix for LLRA to specify different
parameter restrictions in \code{LLRA}.    
}
\usage{
collapse_W(W, listItems, newNames)
}
\arguments{
  \item{W}{A design matrix (for LLRA), typically from a call to
    \code{\link{build_W}} or component \code{$W} from \code{\link{LLRA}}
    or \code{\link{LPCM}}
  } 
  \item{listItems}{A list of numeric vectors. Each component of the list specifies
    columns to be collapsed together. 
  }
  \item{newNames}{An (optional) character vector specifying the names of
    the collapsed effects. 
  }
}
\details{
This function is a convenience function to collapse a design matrix,
i.e. to specify linear trend or treatment effects and so on. Collapsing
here means that effects in columns are summed up. For this, a list of numeric
vectors with the column indices of columns to be collapsed have to be
passed to the function. For example, if you want to collapse column 3, 6
and 8 into one new effect and 1, 4 and 9 into another it needs to be
passed with \code{list(c(3,6,8),c(1,4,9))}.

The new effects can be given names by passing a character vector to the
function with equal length as the list. 
}

\value{An LLRA design matrix as described by Hatzinger and Rusch
  (2009). This can be passed as the \code{W} argument to \code{LLRA} or
 \code{LPCM}.
}
\references{
Hatzinger, R. and Rusch, T. (2009) IRT models with relaxed assumptions
in eRm: A manual-like instruction. \emph{Psychology Science Quarterly},
  \bold{51}, pp. 87--120, \url{http://erm.r-forge.r-project.org/psq_1_2009_06_87-120.pdf}
}
\author{Thomas Rusch}
\seealso{
The function to build design matrices from scratch, \code{\link{build_W}}.
}
\examples{
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.    
llraDat2a <- matrix(unlist(llraDat2[1:20]),ncol=4)
groupvec <-rep(1:3*5,each=20)
W <- build_W(llraDat2a, nitems=5, mpoints=4, grp_n=c(10,20,40), groupvec=groupvec,itmgrps=1:5)

#There are 55 parameters to be estimated
dim(W)

#Imposing a linear trend for the second item ,i.e. parameters in
#columns 32, 37  and 42 need to be collapsed into a single column. 
collItems1 <- list(c(32,37,42))
newNames1 <- c("trend.I2")
Wstar1 <- collapse_W(W,collItems1)

#53 parameters need to be estimated
dim(Wstar1)
}