\encoding{UTF-8}
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
\description{\code{pmat} computes the theoretical person-item matrix with solving
probabilities for each category (except 0th). \code{residuals} computes the squared and standardized residuals based on
the observed and the expected person-item matrix. Chi-square based itemfit and personfit
statistics can be obtained by using \code{itemfit} and \code{personfit}.
}
\usage{
\method{pmat}{ppar}(object)
\method{residuals}{ppar}(object,...)
\method{itemfit}{ppar}(object)
\method{personfit}{ppar}(object)
\method{print}{ifit}(x, visible = TRUE, ...)
\method{print}{pfit}(x, visible = TRUE, ...)
\method{print}{resid}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{ppar}, derived from \code{person.parameter}.}
  \item{x}{Object of class \code{ifit}, \code{pfit}, or \code{resid}.}
  \item{visible}{if \code{FALSE}, returns the matrix of fit statistics that otherwise would be printed.}
  \item{...}{Further arguments passed to or from other methods. They are ignored in this function.}
}
%\details{}
\value{
%Function \code{pmat}:
  \item{pmat}{Matrix of theoretical probabilities for each category except 0th (from function \code{pmat}).}

%Function \code{itemfit} returns a list of class \code{ifit} with components:
  \item{i.fit}{Chi-squared itemfit statistics (from function \code{itemfit}).}
  \item{i.df}{Degrees of freedom for itemfit statistics (from function \code{itemfit}).}
  \item{st.res}{Standardized residuals (from function \code{itemfit}).}
  \item{i.outfitMSQ}{Outfit mean-square statistics (from function \code{itemfit}).}
  \item{i.infitMSQ}{Infit mean-square statistics (from function \code{itemfit}).}

%Function \code{personfit} returns a list of class \code{pfit} with components:
  \item{p.fit}{Chi-squared personfit statistics (from function \code{personfit}).}
  \item{p.df}{Degrees of freedom for personfit statistics (from function \code{personfit}).}
  \item{st.res}{Standardized residuals (from function \code{personfit}).}
  \item{p.outfitMSQ}{Outfit mean-square statistics (from function \code{personfit}).}
  \item{p.infitMSQ}{Infit mean-square statistics (from function \code{personfit}).}

}
\references{
Smith Jr., E. V., and Smith, R. M. (2004). Introduction to Rasch Measurement.
JAM press.

Wright, B.D., and Masters, G.N. Computation of OUTFIT and INFIT Statistics.
Rasch Measurement Transactions, 1990, 3:4 p.84-5
}
\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{person.parameter}}
}
\examples{

# Rasch model, estimation of item and person parameters
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
