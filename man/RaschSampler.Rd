\encoding{UTF-8}
\name{RaschSampler}
\alias{RaschSampler}
%\docType{package}
\title{Rasch Sampler Package}
\description{
The package implements an \acronym{MCMC} algorithm for sampling of
binary matrices with fixed margins complying to the Rasch model.
Its stationary distribution is uniform. The algorithm also allows
for square matrices with fixed diagonal.\cr

Parameter estimates in the Rasch model only depend on the marginal totals of
the data matrix that is used for the estimation. From this it follows that, if the
model is valid, all binary matrices with the same marginals as the observed one
are equally likely. For any statistic of the data matrix, one can approximate
the null distribution, i.e., the distribution if the Rasch model is valid, by taking
a random sample from the collection of equally likely data matrices and constructing
the observed distribution of the statistic.
One can then simply determine the exceedence probability of the statistic in the
observed sample, and thus construct a non-parametric test of the Rasch model.
The main purpose of this package is the implementation of a methodology to build nonparametric
tests for the Rasch model. \cr

In the context of social network theories, where the structure of binary asymmetric
relations is studied, for example,
person \eqn{a}{a} esteems person \eqn{b}{b}, which correponds to a 1 in cell \eqn{(a, b)}{(a, b)}
of the associated adjacency matrix. If one wants to study
the distribution of a statistic defined on the adjacency matrix and conditional
on the marginal totals, one has to exclude the diagonal cells from consideration, i.e.,
by keeping the diagonal cells fixed at an arbitrary value. The \code{RaschSampler} package
has implemented an appropriate option, thus it can be also used for sampling random adjacency
matrices with given marginal totals.
}
\details{
\tabular{ll}{%
Package:\tab RaschSampler\cr
Type:   \tab Package\cr
Version:\tab 0.8-6\cr
Date:   \tab 2012-07-03\cr
License:\tab GNU GPL 2, June 1991\cr%
}
The user has to supply a binary input matrix. After defining appropriate control
parameters using \code{\link{rsctrl}} the sampling function \code{\link{rsampler}}
may be called to obtain an object of class \code{\link{RSmpl}} which contains the
generated random matrices in encoded form. After defining an appropriate function
to operate on a binary matrix (e.g., calculate a statistic such as \code{\link{phi.range}})
the application of this function to the sampled matrices is performed
using \code{\link{rstats}}. Prior to applying the user defined function, \code{\link{rstats}}
decodes the matrices packed in the \code{\link{RSmpl}}-object.\cr

The package also defines a utility function \code{\link{rsextrobj}} for extracting certains parts from
the \code{\link{RSmpl}}-object resulting in an object of class \code{\link{RSmplext}}.
Both types of objects can be saved and reloaded for later use.\cr

Summary methods are available to print information on these objects, as well as
on the control object \code{\link{RSctr}} which is obtained from using
\code{\link{rsctrl}} containing the specification for the sampling routine.

}
\author{Reinhold Hatzinger, Patrick Mair, Norman D. Verhelst}
\references{
Verhelst, N. D. (2008) An Efficient \acronym{MCMC} Algorithm to Sample Binary Matrices with Fixed Marginals. Psychometrika, Volume 73, Number 4\cr
Verhelst, N. D., Hatzinger, R., and Mair, P. (2007) The Rasch Sampler, Journal of Statistical Software, Vol. 20, Issue 4, Feb 2007
}
\note{The current implementation allows for data matrices up to 4096 rows and 128 columns.
      This can be changed by setting \code{nmax} and \code{kmax} in \code{RaschSampler.f90}
      to values which are a power of 2. These values should also be changed in \code{rserror.R}.

      For convenience, we reuse the Fortran code of package version 0.8-1 which cicumvents the
      compiler bug in Linux distributions of GCC 4.3. The following note from package version 0.8-3
      is thus obsolete:
      In case of compilation errors (due to a bug in Linux distributions of GCC 4.3) please use
      \code{RaschSampler.f90} from package version 0.8-1 and change \code{nmax} and \code{kmax}
      accordingly (or use GCC 4.4).}
\keyword{package}
