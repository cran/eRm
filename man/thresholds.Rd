\encoding{UTF-8}
\name{thresholds}
\alias{thresholds}
\alias{thresholds.eRm}
\alias{print.threshold}
\alias{summary.threshold}
\alias{confint.threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computation of item-category treshold parameters.}
\description{This function transforms the beta parameters into threshold
parameters. These can be interpreted by means of log-odds as visualized in ICC plots.
}
\usage{
\method{thresholds}{eRm}(object)
\method{print}{threshold}(x, ...)
\method{summary}{threshold}(object, ...)
\method{confint}{threshold}(object, parm, level = 0.95, ...)
}

\arguments{
Arguments for \code{thresholds}:
  \item{object}{Object of class \code{eRm}.}

Arguments for \code{print}, \code{summary}, and \code{confint} methods:
  \item{x}{Object of class \code{threshold}.}
  \item{parm}{Parameter specification (ignored).}
  \item{level}{Alpha-level.}
  \item{...}{Further arguments to be passed to methods. They are ignored.}
}
\details{For dichotomous models (i.e., RM and LLTM) threshold parameters are not computed.
The \code{print} method returns a location parameter for each item which is the
mean of the corresponding threshold parameters. For LPCM and LRSM the thresholds are
computed for each design matrix block (i.e., measurement point/group) separately
(PCM and RSM have only 1 block).}

\value{
The function \code{thresholds} returns an object of class \code{threshold} containing:
  \item{threshpar}{Vector with threshold parameters.}
  \item{se.thresh}{Vector with standard errors.}
  \item{threshtable}{Data frame with location and threshold parameters.}
}

\references{
Andrich, D. (1978). Application of a psychometric rating model to ordered categories which are scored with successive integers. Applied Psychological Measurement, 2, 581-594.
}
\seealso{
      \code{\link{plotICC.Rm}}
}
\examples{
#Threshold parameterization for a rating scale model
res <- RSM(rsmdat)
th.res <- thresholds(res)
th.res
confint(th.res)
summary(th.res)

#Threshold parameters for a PCM with ICC plot
res <- PCM(pcmdat)
th.res <- thresholds(res)
th.res
plotICC(res)

#Threshold parameters for a LPCM:
#Block 1: t1, g1; Block 2: t1, g2; ...; Block 6: t2,g3
G <- c(rep(1,7),rep(2,7),rep(3,6)) # group vector for 3 groups
res <- LPCM(lpcmdat, mpoints = 2, groupvec = G)
th.res <- thresholds(res)
th.res
}
\keyword{models}
