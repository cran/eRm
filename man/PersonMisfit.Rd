\name{PersonMisfit}
\alias{PersonMisfit}
\alias{PersonMisfit.ppar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Person Misfit}
\description{This function counts the number of persons who do not fit the Rasch model. More specifically, it returns the proportion and frequency of persons - or more generally cases - who exceed a	Chi-square based Z-value of 1.96 (suggesting a statistically significant deviation from the predicted response pattern).
}
\usage{
\method{PersonMisfit}{ppar}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{ppar}.}
}
\details{
  Returns the proportion and absolute number of persons who do not fit the Rasch model (Z-values > 1.96).  
}

\value{
\code{PersonMisfit} returns an object of class \code{MisfittingPersons} containing:
  \item{PersonMisfit}{the proportion of misfitting persons,}
  \item{count.misfit.Z}{the frequency of misfitting person,}
  \item{total.persons}{the number of persons for whom a fit value was estimated.}
}


\author{Adrian Bruegger}
%\note{}


\examples{
rm1 <- RM(raschdat1)
pers <- person.parameter(rm1)
pmfit <- PersonMisfit(pers)
pmfit
summary(pmfit)
}
\keyword{person misfit}

