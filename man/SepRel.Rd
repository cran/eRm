\encoding{UTF-8}
\name{Separation Reliability}
\alias{SepRel}
\alias{print.eRm_SepRel}
\alias{summary.eRm_SepRel}
%
%
%
\title{Person Separation Reliability}
%
\description{%
  This function calculates the proportion of person variance that is not due to error.
  The concept of person separation reliability is very similar to reliability indices such as Cronbach's \eqn{\alpha}.
}
%
%
%
\usage{SepRel(pobject)

\method{print}{eRm_SepRel}(x, \dots)

\method{summary}{eRm_SepRel}(object, \dots)}
%
\arguments{
  \item{pobject}{Object of class \code{ppar} (see \code{\link{person.parameter}}).}
  \item{x}{Object of class \code{eRm_SepRel}.}
  \item{object}{Object of class \code{eRm_SepRel}.}
  \item{\dots}{Further arguments.}
}
%
%
%
\details{
  Returns the person separation reliability \eqn{\frac{\mathrm{SSD}-\mathrm{MSE}}{\mathrm{SSD}}}{(SSD-MSE)/SSD} where SSD is the squared standard deviation and MSE the mean squared error. Note that persons with raw scores of 0 or k are ignored in the computation. 
  \subsection{Caveats}{%
    Please note that the concept of \emph{reliability} and associated problems are fundamentally different between \acronym{IRT} and \acronym{CTT} (Classical Test Theory).
    Separation reliability is more like a workaround to make the \dQuote{change} from \acronym{CTT} to \acronym{IRT} easier for users by providing something \dQuote{familiar.}
    Hence, we recommend not to put too much emphasis on this particular measure and use it with caution.
  }
  \subsection{Varying results in different programs}{%
    If you compare the separation reliability obtained using \pkg{eRm} with values by other software, you will find that they are most likely not equal.
    This has a couple of reasons, one of the most important is the employed estimation method.
    
    \pkg{eRm} uses a conditional maximum likelihood (\acronym{CML}) framework and handles missing values as separate groups during the estimation of item parameters.
    Person parameters are computed in a second step using unconditional or joint maximum likelihood (\acronym{UML} or \acronym{JML}) estimation with item parameters assumed to be known from the first step.
    Other programs might do \acronym{JML} to estimate item and person parameters at the same time, or employ marginal maximum likelihood \acronym{MML} to estimate item parameters, assuming a certain distribution for person parameters.
    In the latter case person parameters might be obtained by various methods like \acronym{EAP}, \acronym{MAP}, \ldots.
    Even \acronym{CML}-based programs yield different values, for example, if they use Warm's weighted maximum likelihood estimation \acronym{WLE} to compute person parameters in the second step.
    
    The bottom line is that, since there is not \dQuote{definite} solution for this problem, you will end up with different values under different circumstances.
    This is another reason to take results and implications with a grain of salt.
  }
}
%
%
%
\value{\code{SepRel} returns a list object of class \code{eRm_SepRel} containing:
  \item{sep.rel}{the person separation reliability,}
  \item{SSD.PS}{the squared standard deviation (i.e., total person variability),}
  \item{MSE}{the mean square measurement error (i.e., model error variance).}
}
%
%
%
\references{%
  Wright, B.D., and Stone, M.H. (1999). \emph{Measurement essentials.} Wide Range Inc., Wilmington. (\url{https://www.rasch.org/measess/me-all.pdf} 28Mb).%
}
\author{Original code by Adrian Brügger (\email{Adrian.Bruegger@imu.unibe.ch}), adapted by Marco J. Maier}
%\note{}
%
%
%
\examples{# Compute Separation Reliability for a Rasch Model:
pers <- person.parameter(RM(raschdat1))
res <- SepRel(pers)
res
summary(res)
}
