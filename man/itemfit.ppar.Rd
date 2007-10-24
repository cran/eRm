\name{itemfit.ppar}
\alias{itemfit.ppar}
\alias{itemfit}
\alias{personfit.ppar}
\alias{personfit}
\alias{residuals.ppar}
\alias{pmat.ppar}
\alias{pmat}
\alias{print.ifit}
\alias{print.pfit}
\alias{print.resid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Residuals, Personfit and Itemfit Statistics}
\description{\code{pmat} computes the expected person-item matrix with corresponding
probabilities. \code{residuals} computes the squared and standardized residuals based on
the observed and the expected person-item matrix. Chi-square based itemfit and personfit
statistics can be obtained by using \code{itemfit} and \code{personfit}.
}
\usage{
\method{pmat}{ppar}(object)
\method{itemfit}{ppar}(object)
\method{personfit}{ppar}(object)
\method{residuals}{ppar}(object,...)
\method{print}{ifit}(x, ...)
\method{print}{pfit}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{ppar}, derived from \code{person.parameter}.}
  \item{x}{Object of class \code{ifit}, \code{pfit}, or \code{resid}.}
  \item{...}{Further arguments passed to or from other methods. They are ignored in this function.}
}
%\details{}
\value{
Function \code{pmat}:
  \item{pmat}{Matrix of theoretical probabilities.}

Function \code{itemfit} returns a list of class \code{ifit} with components:
  \item{i.fit}{Itemfit statistics.}
  \item{i.df}{Degrees of freedom for itemfit statistics.}
  \item{st.res}{Standardized residuals.}

Function \code{personfit} returns a list of class \code{pfit} with components:
  \item{p.fit}{Personfit.}
  \item{p.df}{Degrees of freedom for personfit statistics.}
  \item{st.res}{Standardized residuals.}

}
\references{
Smith Jr., E. V., and Smith, R. M. (2004). Introduction to Rasch Measurement.
JAM press.
}
\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{person.parameter}}
}
\examples{

# Rasch model, estimation of item and person parameters
data(raschdat2)
res <- RM(raschdat2)
p.res <- person.parameter(res)

# Matrix with expected probabilities and corresponding residuals
pmat(p.res)
residuals(p.res)

#Itemfit
itemfit(p.res)

#Personfit
personfit(p.res)

}
\keyword{models}
