\encoding{UTF-8}
\name{llraDat1}
\alias{llraDat1}
\docType{data}
\title{An Artificial LLRA Data Set}
\description{Artificial data set of 5 items, 5 time points and 5 groups for LLRA.}
\usage{llraDat1}
\format{
  A data frame with 150 observations of 26 variables.
  \describe{
    \item{\code{t1.I1}}{Answers to item 1 at time point 1}
    \item{\code{t1.I2}}{Answers to item 2 at time point 1}
    \item{\code{t1.I3}}{Answers to item 3 at time point 1}
    \item{\code{t1.I4}}{Answers to item 4 at time point 1}
    \item{\code{t1.I5}}{Answers to item 5 at time point 1}
    \item{\code{t2.I1}}{Answers to item 1 at time point 2}
    \item{\code{t2.I2}}{Answers to item 2 at time point 2}
    \item{\code{t2.I3}}{Answers to item 3 at time point 2}
    \item{\code{t2.I4}}{Answers to item 4 at time point 2}
    \item{\code{t2.I5}}{Answers to item 5 at time point 2}
    \item{\code{t3.I1}}{Answers to item 1 at time point 3}
    \item{\code{t3.I2}}{Answers to item 2 at time point 3}
    \item{\code{t3.I3}}{Answers to item 3 at time point 3}
    \item{\code{t3.I4}}{Answers to item 4 at time point 3}
    \item{\code{t3.I5}}{Answers to item 5 at time point 3}
    \item{\code{t4.I1}}{Answers to item 1 at time point 4}
    \item{\code{t4.I2}}{Answers to item 2 at time point 4}
    \item{\code{t4.I3}}{Answers to item 3 at time point 4}
    \item{\code{t4.I4}}{Answers to item 4 at time point 4}
    \item{\code{t4.I5}}{Answers to item 5 at time point 4}
    \item{\code{t5.I1}}{Answers to item 1 at time point 5}
    \item{\code{t5.I2}}{Answers to item 2 at time point 5}
    \item{\code{t5.I3}}{Answers to item 3 at time point 5}
    \item{\code{t5.I4}}{Answers to item 4 at time point 5}
    \item{\code{t5.I5}}{Answers to item 5 at time point 5}
    \item{\code{groups}}{The group membership}
  }
}
\details{
This is a data set as described in Hatzinger and Rusch (2009). 5 items
were measured at 5 time points (in columns). Each row corresponds to one
person (P1 to P150).
There are 4 treatment groups and a control group. Treatment group G5 has
size 10 (the first ten subjects),
treatment group G4 has size 20, treatment group G3 has size 30, treatment
group G2 has size 40 and the control group CG has size 50 (the last 50
subjects). Item 1 is dichotomous, all others are polytomous. Item 2, 3, 4 and 5 have 3, 4, 5, 6 categories
respectively.   
}
\references{
Hatzinger, R. and Rusch, T. (2009) IRT models with relaxed assumptions
in eRm: A manual-like instruction. \emph{Psychology Science Quarterly}, \bold{51},
pp. 87--120.
}
\examples{
llraDat1
}
\keyword{datasets}
