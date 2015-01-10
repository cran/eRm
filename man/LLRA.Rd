\encoding{UTF-8}
\name{LLRA}
\alias{LLRA}
\alias{print.llra}

\title{Fit Linear Logistic Models with Relaxed Assumptions (LLRA)
}
\description{
Automatically builds design matrix and fits LLRA.    
}
\usage{
LLRA(X, W, mpoints, groups, baseline, itmgrps = NULL, ...)

\method{print}{llra}(x, ...)
}
\arguments{
  \item{X}{Data matrix as described in Hatzinger and Rusch (2009). It
    must be of wide format, e.g. for each person all item answers are
    written in columns for t1, t2, t3 etc. Hence each row corresponds to
    all observations for a single person. See llraDat1 for an example.
    Missing values are not allowed.    
  }
  \item{W}{Design Matrix for LLRA to be passed to \code{LPCM}. If missing, it
    is generated automatically.
  }
  \item{mpoints}{The number of time points.
  }
  \item{groups}{Vector, matrix or data frame with subject/treatment
    covariates.
  }
  \item{baseline}{An optional vector with the baseline values for the
    columns in group.
  }  
  \item{itmgrps}{
Specifies how many groups of items there are. Currently not functional but may be useful in the future.
}
  \item{x}{For the print method, an object of class \code{"llra"}.
} 
  \item{\dots}{
Additional arguments to be passed to and from other methods. 
}
}
\details{The function \code{LLRA} is a wrapper for \code{LPCM} to fit
  Linear Logistic Models with Relaxed Assumptions (LLRA). LLRA
  are extensions of the LPCM for the measurement of change over a number
  of discrete time points for a set of
  items. It can incorporate categorical covariate information. If no
  design matrix W is passed as an argument, it is built automatically
  from scratch.

  Unless passed by the user, the baseline group is always the one with
  the lowest (alpha-)numerical value for argument \code{groups}. All
  other groups  are labeled decreasingly according to the
  (alpha)-numerical value, e.g. with 2 treatment groups (TG1 and TG2)
  and one control group (CG), CG will be the baseline than TG1 and TG2.
  Hence the group effects are ordered like
  \code{rev((unique(names(groupvec)))} for naming.    

  Caution is advised as LLRA will fail if all changes for a group will be into a
  single direction (e.g. all subjects in the treatment group show
  improvement). Currently only data matrices are supported as arguments.    
}
\value{
Returns an object of class \code{'llra'} (also inheriting from class \code{'eRm'}) containing

\item{loglik}{Conditional log-likelihood.}
\item{iter}{Number of iterations.}
\item{npar}{Number of parameters.}
\item{convergence}{See code output in nlm.}
\item{etapar}{Estimated basic item parameters. These are the LLRA
effect parameters.}
\item{se.eta}{Standard errors of the estimated basic item parameters.}
\item{betapar}{Estimated item (easiness) parameters of the virtual
items (not useful for interpretation here).}
\item{se.beta}{Standard errors of virtual item parameters (not useful for interpretation here).}
\item{hessian}{Hessian matrix if \code{se = TRUE}.}
\item{W}{Design matrix.}
\item{X}{Data matrix in long format. The columns correspond to the
  measurement points and each persons item answers are listed
  susequently in rows.}
\item{X01}{Dichotomized data matrix.}
\item{groupvec}{Assignment vector.}
\item{call}{The matched call.}
\item{itms}{The number of items.}
}
\references{
Fischer, G.H. (1995) Linear logistic models for change. In G.H. Fischer
and I. W. Molenaar (eds.), \emph{Rasch models: Foundations, recent
developments and applications} (pp. 157--181), New York: Springer.

Glueck, J. and Spiel, C. (1997) Item response models for repeated
measures designs: Application and limitations of four different
approaches. \emph{Methods of Psychological
  Research}, \bold{2}. \url{http://www.dgps.de/fachgruppen/methoden/mpr-online/issue2/art6/article.html}

Hatzinger, R. and Rusch, T. (2009) IRT models with relaxed assumptions
in eRm: A manual-like instruction. \emph{Psychology Science Quarterly}, \bold{51},
pp. 87--120, \url{http://erm.r-forge.r-project.org/psq_1_2009_06_87-120.pdf}
}
\author{Thomas Rusch}
\section{Warning}{A warning is printed that the first two categories
  for polytomous items are equated to save parameters. See Hatzinger and
  Rusch (2009) for a justification why this is valid also from a substantive
  point of view.} 
\seealso{
The function to build the design matrix \code{\link{build_W}}, and the
S3 methods \code{\link{summary.llra}} and \code{\link{plotTR}} and
\code{\link{plotGR}} for plotting. 
}
\examples{
##Example 6 from Hatzinger & Rusch (2009)
groups <- c(rep("TG",30),rep("CG",30))
llra1 <- LLRA(llradat3,mpoints=2,groups=groups)
llra1

\dontrun{
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.
dats <- llraDat2[1:20]
groups <- llraDat2$group
tps <- 4

#baseline CG
ex2 <- LLRA(dats,mpoints=tps,groups=groups) 

#baseline TG1
ex2a <- LLRA(dats,mpoints=tps,groups=groups,baseline="TG1") 

#summarize results
summary(ex2)
summary(ex2a)

#plotting
plotGR(ex2)
plotTR(ex2)}
}
