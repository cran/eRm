\encoding{UTF-8}
\name{eRm.data}
\alias{eRm.data}
\alias{raschdat1}
\alias{raschdat1_RM_fitted}
\alias{raschdat1_RM_plotDIF}
\alias{raschdat1_RM_lrres2}
\alias{raschdat2}
\alias{raschdat3}
\alias{raschdat4}
\alias{lltmdat1}
\alias{lltmdat2}
\alias{pcmdat}
\alias{pcmdat2}
\alias{lpcmdat}
\alias{rsmdat}
\alias{lrsmdat}
\docType{data}
\title{Data for Computing Extended Rasch Models}
\description{Artificial data sets for computing extended Rasch models.}
\usage{raschdat1
raschdat2
raschdat3
raschdat4

lltmdat1
lltmdat2

rsmdat

lrsmdat

pcmdat
pcmdat2

lpcmdat

raschdat1_RM_fitted
raschdat1_RM_plotDIF
raschdat1_RM_lrres2}
\format{Numeric matrices with subjects as rows, items as columns, missing values as \code{NA}.}
\details{\code{raschdat1_RM_fitted} is the resulting object of \code{RM(raschdat1)} and used in examples to reduce computation time. For the generation of \code{raschdat1_RM_plotDIF} see the excluded example code of \code{\link{plotDIF}}. \code{raschdat1_RM_lrres2} results from \code{LRtest(RM(raschdat1), split = "mean")}}
\keyword{datasets}
