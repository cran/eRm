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
and its linear extension (LPCM). The parameters are estimated by conditional maximum 
likelihood (CML). Missing values are allowed in the data matrix. Additional features 
are the estimation of the person parameters, LR-Model test, item-spefific Wald test, 
itemfit and personfit statistics, various ICC plots. An eRm platform is provided at 
http://r-forge.r-project.org.
}
\details{
\tabular{ll}{
Package: \tab eRm\cr
Type: \tab Package\cr
Version: \tab 0.9.1.\cr
Date: \tab 2007-05-17\cr
License: \tab GPL\cr
}
The basic input units for the functions are the person-item matrix X and the design matrix W. 
Missing values in X are coded with \code{NA}. 
By default, W is generated automatically, but it can be specified by the user as well. 
The function call of the basic models can be achieved through \code{RM(X, W)}, 
\code{RSM(X, W)}, and \code{PCM(X, W)}.

The linear extensions provide the possibility to fit a more restricted model than its basic complement, but
also a generalization by imposing repeated measurement designs and group contrasts. These models can 
be estimated by using, e.g., the call \code{LLTM(X, W, mpoints = 2, groupvec = G)}, \code{LRSM(X, W, mpoints = 2, groupvec = G)}, and
\code{LPCM(X, W, mpoints = 2, groupvec = G)}. \code{G} is a vector with the group membership for each subject
ordered according to the rows of the data matrix.

\code{RM} produces an object belonging to the classes \code{dRM}, \code{RM}, and 
\code{eRm}. \code{PCM} and \code{RSM} produce objects belonging to the classes 
\code{RM} and \code{eRm}, whereas results of \code{LLTM}, \code{LRSM}, and \code{LLTM} are object of class \code{eRm}.
 
}
\author{Patrick Mair, Reinhold Hatzinger

Maintainer: Patrick Mair <patrick.mair@wu-wien.ac.at>
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for 
the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models 
with the eRm package in R. Psychology Science, 49, 26-43.
}
\keyword{models}

