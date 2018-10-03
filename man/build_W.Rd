\encoding{UTF-8}
\name{build_W}
\alias{build_W}
\alias{build_catdes}
\alias{build_trdes}
\alias{build_effdes}
\alias{get_item_cats}

\title{
Automatized Construction of LLRA Design Matrix 
}
\description{
Builds a design matrix for LLRA from scratch.
}
\usage{
build_W(X, nitems, mpoints, grp_n, groupvec, itmgrps)
}
\arguments{
  \item{X}{Data matrix as described in Hatzinger and Rusch (2009). It
    must be of long format, e.g. for each person all item answers are written in subsequent rows. The columns correspond to time
    points. Missing values are not allowed. It can easily be
    constructed from data in wide format with
    \code{matrix(unlist(data),ncol=mpoints)} or from \code{\link{llra.datprep}}.
  } 
  \item{nitems}{The number of items.
}
  \item{mpoints}{The number of time points.
}
 \item{grp_n}{A vector of number of subjects per g+1 groups (e.g. g
    treatment or covariate groups and 1 control or baseline group.
    The sizes must be ordered like the corresponding groups.
  }
 \item{groupvec}{Assignment vector, i.e. which person belongs to which
   treatment/item group 
  }  
  \item{itmgrps}{Specifies how many groups of items there are. 
 }
}
\details{
The function is designed to be modular and calls four internal function
\code{build_effdes} (for treatment/covariate effects), \code{build_trdes} (for trend
effects), \code{build_catdes} (for category parameter design matrix) and
\code{get_item_cats} (checks how many categories each item has). Those functions are not intended to be used by the user.

Labeling of effects also happens in the internal functions. 
}
\value{An LLRA design matrix as described by Hatzinger and Rusch
  (2009). This can be passed as the \code{W} argument to \code{LLRA} or
 \code{LPCM}.

  The design matrix specifies every item to lie on its own
  dimension. Hence at every time point > 1, there are effects for
  each treatment or covariate group as well as trend effects for every
  item. Therefore overall there are items x (groups-1) x (time points-1)
  covariate effect parameters and items x (time points-1) trend
  parameters specified. For polytomous items there also are parameters
  for each category with the first and second category being equated for each item. They
  need not be equidistant. The number of parameters therefore increase
  quite rapidly for any additional time point, item or covariate group.       
}
\references{
Hatzinger, R. and Rusch, T. (2009) IRT models with relaxed assumptions
in eRm: A manual-like instruction. \emph{Psychology Science Quarterly},
  \bold{51}, pp. 87--120.
}
\author{Thomas Rusch}
\section{Warning }{A warning is printed that the first two categories
  for polytomous items are equated.} 

\seealso{
This function is used for automatic generation of the design matrix in \code{\link{LLRA}}.
}
\examples{
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.
llraDat2a <- matrix(unlist(llraDat2[1:20]),ncol=4)
groupvec <-rep(1:3*5,each=20)
W <- build_W(llraDat2a,nitems=5,mpoints=4,grp_n=c(10,20,40),groupvec=groupvec,itmgrps=1:5)

#There are 55 parameters
dim(W)

\dontrun{
#Estimating LLRA by specifiying W
ex2W <- LLRA(llraDat2[1:20],W=W,mpoints=4,groups=llraDat2[21])}
}