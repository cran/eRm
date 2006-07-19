\name{RSM}
\alias{RSM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of rating scale models}
\description{
  This function computes the parameter estimates of a rating scale model for polytomous 
  item responses by using CML estimation. }
\usage{
RSM(X, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix or data frame with item responses (starting from 0); rows represent individuals, columns represent items.}
  \item{W}{Design matrix for the RSM. If omitted, the function will compute W automatically.}
}
\details{
  The design matrix approach transforms the RSM into a partial credit model 
  and estimates the corresponding basic parameters by using CML.  
  Available methods for RSM-objects are \code{print}, \code{coef}, \code{model.matrix}, 
  \code{vcov}, \code{plot}, \code{summary}.
}
\value{
  Returns an object of class \code{Rm} and contains the log-likelihood value, 
  the parameter estimates and their standard errors.
  
  \item{loglik}{The log-likelihood.}
  \item{iter}{Number of iterations required.}
  \item{etapar}{Estimated basic item parameters.}
  \item{se_eta}{Standard errors of the estimated basic item parameters.}
  \item{betapar}{Estimated item parameters.}
  \item{LR}{The log-likelihood test statistic for the model.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{Missing Values are not allowed in X. All items must have the same number of 
response categories. 
}

\seealso{\code{\link{print.eRm}},\code{\link{coef.eRm}},\code{\link{vcov.eRm}},\code{\link{model.matrix.eRm}},\code{\link{plot.Rm}},
\code{\link{summary.Rm}}
}
\examples{

#PCM with 10 subjects, 3 items

data(X_rsm)
res <- RSM(X_rsm)
res              
}

\keyword{models}
