\name{eRm-package}
\alias{eRm-package}
\alias{eRm}
\docType{package}
\title{
extended Rasch modeling
}
\description{
This package estimates extended Rasch models, i.e. the 
ordinary Rasch model for dichotomous data (RM), the linear logistic test model (LLTM),
the rating scale model (RSM) and its linear extension (LRSM), the partial credit model (PCM)
and its linear extension (LPCM). The parameters are estimated by conditional maximum likelihood (CML). 
}
\details{
\tabular{ll}{
Package: \tab eRm\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2006-10-06\cr
License: \tab GPL2 or higher\cr
}
The basic input units for the functions are the person-item-matrix X and the design matrix W. 
By default, W is generated automatically, but it can be specified by the user as well. The function
call of the basic models is \code{RM(X, W)}, \code{RSM(X, W)}, and \code{PCM(X, W)}.

The linear extensions provide the possibility to fit a more restricted model than its basic complement, but
also a generalization by imposing repeated measurement designs and group contrasts. These models can 
be estimated by using e.g. the call \code{LLTM(X, W, mpoints = 2, groupvec = G)}, \code{LRSM(X, W, mpoints = 2, groupvec = G)}, and
\code{LPCM(X, W, mpoints = 2, groupvec = G)}. \code{G} is a vector with the group membership of each subject.
 
}
\author{Patrick Mair, Reinhold Hatzinger

Maintainer: Patrick Mair <patrick.mair@wu-wien.ac.at>
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.
}
\keyword{models}

