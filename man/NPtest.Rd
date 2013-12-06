\name{NPtest}
\Rdversion{1.1}
\alias{NPtest}
\title{function to perform nonparametric Rasch model tests}
\description{A variety of nonparametric tests as proposed by Ponocny(2001), Koller and Hatzinger(2012),
    and an 'exact' version of the Martin-Loef test are implemented. The function operates on
    random binary matrices that have been generated using an
    MCMC algorithm (Verhelst, 2008) from the RaschSampler package (Hatzinger, Mair, and Verhelst, 2009).
}
\usage{
NPtest(obj, n = NULL, method = "T1", ...)
}
\arguments{
  \item{obj}{
     A binary data matrix (or data frame) or
     an object containing the output from the \code{\link{RaschSampler}} package.
  }
  \item{n}{
     If \code{obj} is a matrix or a data frame, \code{n} n is the number of sampled matrices
     (default is 500)
  }
  \item{method}{
     One of the test statistics. See details below.
  }
  \item{\dots}{
     Further arguments according to \code{method}. See details below. Additionally, the sampling routine
     can be controlled by specifying \code{burn_in =}, \code{step =}, and \code{seed =} (for details see
     below and . A summary of the sampling object may be obtained using
     the option \code{RSinfo = TRUE}.
  }
}
\details{
     The function uses the  \code{\link{RaschSampler}} package. It has to be installed to use
     \code{NPtest}.
     On input the user has to supply
     either a binary data matrix or a RaschSampler output object. If the input is a data matrix, the RaschSampler
     is called with default values (i.e., \code{rsctrl(burn_in = 256, n_eff = n, step = 32)}, see \code{\link{rsctrl}}),
     where \code{n} corresponds to \code{n_eff} (the default number of sampled matrices is 500).
     By default, the starting values for the random number generators (\code{seed}) are chosen randomly using system time.
     Methods other than those listed below can easily be implemented using the RaschSampler package directly.

     The currently implemented methods (following Ponocny's notation of \emph{T}-statistics) and their options are:
   \describe{
     \item{\bold{T1:}}{\code{method = "T1"}, no further option}\cr
     Checks for local dependence via increased inter-item correlations. For all item pairs
     cases are counted with equal responses on both items.

     \item{\bold{T1m:}}{\code{method = "T1m"}, no further option}\cr
     Checks for multidimensionality via decreased inter-item correlations. For all item pairs
     cases are counted with equal responses on both items.

     \item{\bold{T1l:}}{\code{method = "T1l"}, no further option}\cr
     Checks for learning. For all item pairs
     cases are counted with response pattern (1,1).

     \item{\bold{Tmd:}}{\code{method = "Tmd", idx1, idx2}}\cr
     \code{idx1, idx2} \ldots vectors of indices specifying items which define two subscales, e.g.,
     \code{idx1 = c(1, 5, 7)} and \code{idx2 = c(3, 4, 6)} \cr
     Checks for multidimensionality based on correlations of person raw scores for the subscales.

     \item{\bold{T2:}}{\code{method = "T2", idx = NULL, stat = "var"}}\cr
     \code{idx} \ldots vector of indices specifying items which define a subscale, e.g., \code{idx = c(1, 5, 7)}\cr
     \code{stat} \ldots one of \code{"var"} (variance), \code{"mad1"} (mean absolute deviation),
     \code{"mad2"} (median absolute deviation), \code{"range"} (range)\cr
     Checks for local dependence within model deviating subscales via increased
     dispersion of subscale person rawscores.

     \item{\bold{T2m:}}{\code{method = "T2m", idx = NULL, stat = "var"}}\cr
     \code{idx} \ldots vector of indices specifying items which define a subscale, e.g., \code{idx = c(1, 5, 7)}\cr
     \code{stat} \ldots one of \code{"var"} (variance), \code{"mad1"} (mean absolute deviation),
     \code{"mad2"} (median absolute deviation), \code{"range"} (range)\cr
     Checks for multidimensionality within model deviating subscales via decreased
     dispersion of subscale person rawscores.

     \item{\bold{T4:}}{\code{method = "T4", idx = NULL, group = NULL, alternative = "high"}}\cr
     \code{idx} \ldots vector of indices specifying items which define a subscale, e.g., \code{idx = c(1, 5, 7)}\cr
     \code{group} \ldots logical vector defining a subject group, e.g., \code{group = (age >= 15 && age < 30)}\cr
     \code{alternative} \ldots one of \code{"high"} or \code{"low"}. Specifies the alternative hypothesis.\cr
     Checks for group anomalies (DIF) via too high (low) raw scores on item(s) for specified group.

%  removed in version 0.14-5
%     \item{\bold{T7:}}{\code{method = "T7", idx = NULL}}\cr
%     \code{idx} \ldots vector of indices specifying items which define a subscale, e.g., \code{idx = c(1, 5, 7)}\cr
%     Checks for lower discrimination (2PL) in item subscale via counting cases with response 1 on more
%     difficult and 0 on easier items. The test is global for the subscale, i.e. all subscale items are evaluated
%     using a single statistic.
%
%     \item{\bold{T7a:}}{\code{method = "T7a", idx = NULL}}\cr
%     \code{idx} \ldots vector of indices specifying items to investigate, e.g., \code{idx = c(1, 5, 7)}\cr
%     Checks for lower discrimination (2PL) of an item compared to another item  via counting cases with response 1 on more
%     difficult and 0 on easier item. The test is performed pairwise, i.e. a statistic is calculated for each item pair.
%     # Examples
%     ##---- T7, T7a --------------------------------------------------
%     # simultaenous test for all items in subscale
%     t7<-NPtest(rmat,method="T7",idx=1:3)
%     t7
%
%     # test for item-pairs
%     t7a<-NPtest(rmat,method="T7a",idx=c(1,3,5)) # test for item-pairs
%     t7a

     \item{\bold{T10:}}{\code{method = "T10", splitcr="median"}}\cr
     \code{splitcr} \ldots split criterion for subject raw score splitting. \code{"median"} uses the median as split criterion,
     \code{"mean"} performs a mean-split. Optionally \code{splitcr} can also be a vector which assigns each person to a
     one of two subgroups (e.g., following an external criterion). This vector can be numeric, character, logical or a factor.\cr
     Gobal test for subgroup-invariance. Checks for different item difficulties in two subgroups (for details see Ponocny, 2001).

     \item{\bold{T11:}}{\code{method = "T11"}, no further specific option}\cr
     Gobal test for local dependence. The statistic calculates the sum of absolute deviations between the observed inter-item correlations
     and the expected correlations.

     \item{\bold{Tpbis:}}{\code{method = "Tpbis", idxt, idxs}}\cr
     Test for discrimination. The statistic calculates a pointbiserial correlation for a test item (specified via \code{idxt})
     with the person row scores for a subscale of the test sum (specified via \code{idxs}). If correlation is too low, the test item
     shows different discrimination compared to the items of the subscale.
   }

    The 'exact' version of the \bold{Martin-Loef} statistic is specified via \code{method = "MLoef"} and optionally \code{splitcr}
    (see \code{\link{MLoef}}).
}
\value{
Depends on the method used. For each method a list is returned. The returned objects are of class
\code{T1obj},\code{T1mobj},\code{T1lobj}, \code{Tmdobj}, \code{T2obj}, \code{T2mobj}, \code{T4obj},
\code{T10obj}, \code{T11obj}, and \code{Tpbisobj} corresponding to the method used.
The main output element is \code{prop} giving the one-sided p-value, i.e., the number of statistics from the sampled matrices which are equal
or exceed the statistic based on the observed data. For \emph{T1}, \emph{T1m}, and \emph{T1l} \code{prop} is a vector.
For the \emph{Martin-Loef} test the returned object is of class \code{MLobj}. Besides other elements, it contains a \code{prop} vector and \code{MLres}, the output
object from the asymptotic Martin-Loef test on the input data.
}
\note{
The \code{RaschSampler} package is no longer required to use \code{NPtest} since version 0.15-0.
}
\references{
Ponocny, I. (2001) Nonparametric goodness-of-fit tests for the rasch model. Psychometrika,  Volume 66, Number 3\cr
Verhelst, N. D. (2008) An Efficient MCMC Algorithm to Sample Binary
Matrices with Fixed Marginals. Psychometrika, Volume 73, Number 4\cr
Verhelst, N. D., Hatzinger, R., and Mair, P. (2007) The Rasch Sampler, Journal of Statistical Software, Vol. 20, Issue 4, Feb 2007
}
\author{
Reinhold Hatzinger
}
%\note{
%Maybe notes appear here
%}
%\seealso{
%  \code{\link{RaschSampler}}
%}
\examples{
### Preparation:

# data for examples below
X <- raschdat1

# generate 100 random matrices based on original data matrix
rmat <- rsampler(X, rsctrl(burn_in = 100, n_eff = 100, seed = 123))

## the following examples can also directly be used by setting
## rmat <- raschdat1
## without calling rsampler() first, e.g.,
t1 <- NPtest(raschdat1, n = 100, method = "T1")


### Examples ###################################################################

###--- T1 ----------------------------------------------------------------------
t1 <- NPtest(rmat, method = "T1")
# choose a different alpha for selecting displayed values
print(t1, alpha = 0.01)


###--- T2 ----------------------------------------------------------------------
t21 <- NPtest(rmat, method = "T2", idx = 1:5, burn_in = 100, step = 20,
              seed = 7654321, RSinfo = TRUE)
# default stat is variance
t21

t22 <- NPtest(rmat, method = "T2", stat = "mad1",
              idx = c(1, 22, 5, 27, 6, 9, 11))
t22


###--- T4 ----------------------------------------------------------------------
age <- sample(20:90, 100, replace = TRUE)
# group MUST be a logical vector
# (value of TRUE is used for group selection)
age <- age < 30
t41 <- NPtest(rmat, method = "T4", idx = 1:3, group = age)
t41

sex <- gl(2, 50)
# group can also be a logical expression (generating a vector)
t42 <- NPtest(rmat, method = "T4", idx = c(1, 4, 5, 6), group = sex == 1)
t42


###--- T10 ---------------------------------------------------------------------
t101 <- NPtest(rmat, method = "T10")       # default split criterion is "median"
t101

\dontrun{
split <- runif(100)
t102 <- NPtest(rmat, method = "T10", splitcr = split > 0.5)
t102

t103 <- NPtest(rmat, method = "T10", splitcr = sex)
t103}


###--- T11 ---------------------------------------------------------------------
t11 <- NPtest(rmat, method = "T11")
t11


###--- Tpbis -------------------------------------------------------------------
tpb <- NPtest(X[, 1:5], method = "Tpbis", idxt = 1, idxs = 2:5)
tpb


###--- Martin-Loef -------------------------------------------------------------
\dontrun{
# takes a while ...
split <- rep(1:3, each = 10)
NPtest(raschdat1, n = 100, method = "MLoef", splitcr = split)}
}
\keyword{htest}
\keyword{nonparametric}
