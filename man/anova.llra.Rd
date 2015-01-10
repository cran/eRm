\encoding{UTF-8}
\name{anova.llra}
\alias{anova.llra}

\title{Analysis of Deviance for Linear Logistic Models with Relaxed Assumptions}
\description{Compute an analysis of deviance table for one or more LLRA.}
\usage{
\method{anova}{llra}(object, ...)
}
\arguments{
  \item{object, ... }{Objects of class "llra", typically the result of a
    call to \code{\link{LLRA}}.
    }
}
\details{
An analysis of deviance table will be calculated. The models in rows are
ordered from the smallest to the largest model. Each row shows the
number of parameters (Npar) and the log-likelihood (logLik). For all but
the first model, the parameter difference (df) and the difference in
deviance or the likelihood ratio (-2LR) is given between two subsequent
models (with increasing complexity). Please note that interpreting these
values only makes sense if the models are nested. 

The table also contains p-values comparing the reduction in the
deviance to the df for each row based on the asymptotic Chi^2-Distribution of the Likelihood ratio test statistic.   
}
\value{
An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\author{Thomas Rusch}
\section{Warning:}{
The comparison between two or more models by \code{anova} will only be valid
if they are fitted to the same dataset and if the models are nested. The
function does not check if that is the case. 
}
\seealso{
The model fitting function \code{\link{LLRA}}.
}
\examples{\dontrun{
##An LLRA with 2 treatment groups and 1 baseline group, 5 items and 4
##time points. Item 1 is dichotomous, all others have 3, 4, 5, 6
##categories respectively.

#fit LLRA
ex2 <- LLRA(llraDat2[,1:20],mpoints=4,groups=llraDat2[,21])

#Imposing a linear trend for items 2 and 3 using collapse_W 
collItems2 <- list(c(32,37,42),c(33,38,43))
newNames2 <- c("trend.I2","trend.I3")
Wnew <- collapse_W(ex2$W,collItems2,newNames2)

#Estimating LLRA with the linear trend for item 2 and 3
ex2new <- LLRA(llraDat2[1:20],W=Wnew,mpoints=4,groups=llraDat2[21])

#comparing models with likelihood ratio test
anova(ex2,ex2new)}}
