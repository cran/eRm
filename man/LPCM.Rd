\name{LPCM}
\alias{LPCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of linear partial credit models}
\description{
  This function computes the parameter estimates of a linear partial credit model (LRSM) 
  for polytomuous item responses by using CML estimation. 
}
\usage{
LPCM(X, W , mpoints = 1, groupvec = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix or data frame with item responses (starting from 0); rows represent individuals (N in total), columns represent items.}
  \item{W}{Design matrix for the LPCM. If omitted, the function will compute W automatically.}
  \item{mpoints}{Number of measurement points.}
  \item{groupvec}{Vector of length N which determines the group membership of each subject, starting from 1}
}
\details{
  Through appropriate definition of W the LPCM can be viewed as a more parsimonous 
  PCM, on the one hand, e.g. by imposing some cognitive base operations 
  to solve the items. One the other hand, linear extensions of the Rasch model 
  such as group comparisons and repeated measurement designs can be computed. 
  If more than one measurement point is examined, the item responses for the 2nd, 3rd, etc. 
  measurement point are added column-wise in X, i.e. X(T1)|X(T2)|... Available methods for LPCM-objects are \code{print}, \code{coef}, 
  \code{model.matrix}, \code{vcov}, \code{summary}.
}
\value{
  Returns on object of class \code{eRm} and contains the log-likelihood value, 
  the parameter estimates and their standard errors.
  
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
  \item{mpoints}{Number of measurement points.}
  \item{ngroups}{Number of groups.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, 
Recent Developements, and Applications. Springer.
}
\author{Patrick Mair, Reinhold Hatzinger}
\note{NA's are not allowed in X, the category coding must start with 0 (lowest category).
}
\seealso{\code{\link{print.eRm}},\code{\link{coef.eRm}},\code{\link{vcov.eRm}},\code{\link{model.matrix.eRm}},\code{\link{summary.eRm}}}
\examples{

#LPCM for two measurement points and two subject groups
#20 subjects, 2*3 items

data(lpcmdat)
G <- c(rep(1,10),rep(2,10))             
res <- LPCM(lpcmdat, mpoints = 2, groupvec = G)
res
}

\keyword{models}
