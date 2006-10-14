\name{RM}
\alias{RM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of Rasch Models}
\description{
  This function computes the parameter estimates of a Rasch model for binary item responses by using CML estimation. 
}
\usage{
RM(X, W)
}
\arguments{
  \item{X}{Input 0/1 data matrix or data frame; rows represent individuals, columns represent items.}
  \item{W}{Design matrix for the Rasch model. If omitted, the function will compute W automatically.}
}
\details{
  For estimating the item parameters the CML method is used. 
  Available methods for RM-objects are \code{print}, \code{coef}, \code{model.matrix}, \code{vcov}, \code{plot}, \code{summary}.
}
\value{
  Returns an object of class \code{Rm} and contains the log-likelihood value, the parameter estimates and their standard errors.
  
  \item{model}{Type of model.}
  \item{loglik}{The log-likelihood.}
  \item{df}{Degrees of freedom.}
  \item{iter}{Number of iterations required.}
  \item{etapar}{Estimated basic item parameters.}
  \item{se_eta}{Standard errors of the estimated basic item parameters.}
  \item{hessian}{Hessian matrix.}
  \item{betapar}{Estimated item parameters.}
  \item{LR}{The log-likelihood test statistic for the model.}
  \item{W}{Design matrix.}
  \item{etaparG1}{Parameters for first LR-group.}
  \item{etaparG2}{Parameters for second LR-group.}

}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{Missing Values are not allowed in X, the category coding must start with 0 (lowest category). 
}
\seealso{\code{\link{print.eRm}},\code{\link{coef.eRm}},\code{\link{vcov.eRm}},\code{\link{model.matrix.eRm}},\code{\link{plot.Rm}},
\code{\link{summary.eRm}}}
\examples{

#Rasch model with 30 items, 100 persons

data(raschdat)
res <- RM(raschdat)
res
summary(res)                
plot(res)                   
}

\keyword{models}
