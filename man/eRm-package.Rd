\encoding{UTF-8}
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
Martin-Loef test, nonparametric Monte-Carlo tests,
itemfit and personfit statistics, various ICC plots.

An \pkg{eRm} platform is provided at \url{http://r-forge.r-project.org/projects/erm/}.
}
\details{%
\tabular{ll}{%
Package: \tab{}eRm\cr%
Type:    \tab{}Package\cr%
Version: \tab{}\Sexpr[stage=build]{packageDescription("eRm")$Version}\cr%
Date:    \tab{}\Sexpr[stage=build]{packageDescription("eRm")$Date}\cr%
License: \tab{}\Sexpr[stage=build]{packageDescription("eRm")$License}\cr%
}
The basic input units for the functions are the person-item matrix X and the design matrix W.
Missing values in X are coded with \code{NA}.
By default, W is generated automatically, but it can be specified by the user as well.
The function call of the basic models can be achieved through \code{RM(X, W)},
\code{RSM(X, W)}, and \code{PCM(X, W)}.

The linear extensions provide the possibility to fit a more restricted model than its basic complement,
such as \code{LLTM(X, W)}, \code{LRSM(X, W)},\code{LPCM(X, W)}, but
also a generalization by imposing repeated measurement designs and group contrasts. These models can
be estimated by using, e.g.,
\code{LLTM(X, W, mpoints = 2, groupvec = g)},\cr
\code{LRSM(X, W, mpoints = 2, groupvec = g)},\cr
\code{LPCM(X, W, mpoints = 2, groupvec = g)},\cr
and as very flexible multidimensional model for repeated measurements
\code{LLRA(X, W, mpoints = 2, groups = G)},\cr
\code{mpoints} specifies the number of measurement or time points,
\code{g} is a vector with the group membership for each subject,
ordered according to the rows of the data matrix, and
\code{G} is a matrix with subject covariates (e.g., treatments),

\code{RM} produces an object belonging to the classes \code{dRm}, \code{Rm}, and
\code{eRm}. \code{PCM} and \code{RSM} produce objects belonging to the classes
\code{Rm} and \code{eRm}, whereas results of \code{LLTM}, \code{LRSM}, \code{LLTM}
and \code{LLRA} are objects of class \code{eRm}. For a detailled overview of all
classes defined in the package and the functions depending on them see the package's vignette.

We acknowledge Julian Gilbey for writing the \code{plotPWmap} function, Kathrin Gruber
for the function \code{plotDIF}, and Thomas Rusch for \code{LLRA}, related utilities and
functionality to calculate and plot item and test information.
The \pkg{eRm} package contains functions from the packages \pkg{sna}, \pkg{gtools} and \pkg{ROCR}.
Thanks to Carter T. Butts, Gregory R. Warnes, and Tobias Sing et al.
}
\note{The fitting engine by default is \code{\link{nlm}} unless changed to \code{\link{optim}}.
      For specification of the optimizer the global variable \code{fitctrl} has to be used, i.e.,
      \code{fitctrl <- "nlm"} or \code{fitctrl <- "optim"}.}


\author{Patrick Mair, Reinhold Hatzinger, Marco J. Maier, and others.

Maintainer: Patrick Mair (\email{mair@fas.harvard.edu})%
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations,
Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The \pkg{eRm} package for
the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models
with the \pkg{eRm} package in R. Psychology Science, 49, 26-43.
}
\keyword{models}
