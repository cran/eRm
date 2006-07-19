\name{LLTM}
\alias{LLTM}
\title{Estimation of linear logistic test models}
\description{
  This function computes the parameter estimates of a linear logistic test model (LLTM) 
  for binary item responses by using CML estimation. 
}
\usage{
LLTM(X, W, mpoints = 1, Groups = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input 0/1 data matrix or data frame; rows represent individuals (N in total), columns represent items.}
  \item{W}{Design matrix for the Rasch model. If omitted, the function will compute W automatically.}
  \item{mpoints}{Number of measurement points.}
  \item{Groups}{Vector of length N which determines the group membership of each subject, starting from 1}
}
\details{
  Through appropriate definition of W the LLTM can be viewed as a more parsimonous 
  Rasch model, on the one hand, e.g. by imposing some cognitive base operations 
  to solve the items. One the other hand, linear extensions of the Rasch model 
  such as group comparisons and repeated measurement designs can be computed. 
  If more than one measurement point is examined, the item responses for the 2nd, 3rd, etc. 
  measurement point are added column-wise in X.   
  Available methods for LLTM-objects are \code{print}, \code{coef}, 
  \code{model.matrix}, \code{vcov}.
}
}
\value{
  Returns on object of class \code{eRm} and contains the log-likelihood value, 
  the parameter estimates and their standard errors.
  
  \item{loglik}{The log-likelihood.}
  \item{iter}{Number of iterations required.}
  \item{etapar}{Estimated basic item parameters.}
  \item{se_eta}{Standard errors of the estimated basic item parameters.}
  \item{betapar}{Estimated item parameters.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{NA's are not allowed in X.
}
\seealso{\code{\link{print.eRm}},\code{\link{coef.eRm}},\code{\link{vcov.eRm}},\code{\link{model.matrix.eRm}}}
\examples{

#LLTM for 2 measurement points 
#100 persons, 2*10 items, W generated automatically

data(X_lltm)                                         
res <- LLTM(X_lltm, mpoints = 2)
res
}

\keyword{models}
