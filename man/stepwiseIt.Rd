\encoding{UTF-8}
\name{stepwiseIt}
\alias{stepwiseIt}
\alias{stepwiseIt.eRm}
\alias{print.step}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepwise item elimination}
\description{This function eliminates items stepwise according to one of the following 
criteria: itemfit, Wald test, Andersen's LR-test
}
\usage{
\method{stepwiseIt}{eRm}(object, criterion = list("itemfit"), alpha = 0.05,
           verbose = TRUE, maxstep = NA)
}

\arguments{
  \item{object}{Object of class \code{eRm}.}
  \item{criterion}{List with either \code{"itemfit"}, \code{"Waldtest"} or \code{"LRtest"} as first element. 
  Optionally, for the Waldtest and LRtest a second element containing the split criterion can be specified (see details).}
  \item{alpha}{Significance level.}
  \item{verbose}{If \code{TRUE} intermediate results are printed out. }
  \item{maxstep}{Maximum number of elimination steps. If \code{NA} the procedure stops when the itemset is Rasch homogeneous.}
}

\details{If \code{criterion = list("itemfit")} the elimination stops when none of the p-values 
in itemfit is significant. Within each step the item with the largest chi-squared 
itemfit value is excluded.

If \code{criterion = list("Waldtest")} the elimination stops when none of the p-values 
resulting from the Wald test is significant. Within each step the item with the largest z-value in 
Wald test is excluded. 

If \code{criterion = list("LRtest")} the elimination stops when Andersen's LR-test is not
significant. Within each step the item with the largest z-value in Wald test is excluded. 
}

\value{
The function returns an object of class \code{step} containing:
  \item{X}{Reduced data matrix (bad items eliminated)}
  \item{fit}{Object of class \code{eRm} with the final item parameter elimination}
  \item{it.elim}{Vector contaning the names of the eliminated items}
  \item{res.wald}{Elimination results for Wald test criterion}
  \item{res.itemfit}{Elimination results for itemfit criterion}
  \item{res.LR}{Elimination results for LR-test criterion}
  \item{nsteps}{Number of elimination steps}
}

\seealso{ \code{\link{LRtest.Rm}}, \code{\link{Waldtest.Rm}}, \code{\link{itemfit.ppar}}
}
\examples{

## 2pl-data, 100 persons, 10 items
set.seed(123)
X <- sim.2pl(500, 10, 0.4)
res <- RM(X)

## elimination according to itemfit
stepwiseIt(res, criterion = list("itemfit"))      

## Wald test based on mean splitting
stepwiseIt(res, criterion = list("Waldtest","mean")) 

## Andersen LR-test based on random split
set.seed(123)
groupvec <- sample(1:3, 500, replace = TRUE)
stepwiseIt(res, criterion = list("LRtest",groupvec))

}
\keyword{models}
