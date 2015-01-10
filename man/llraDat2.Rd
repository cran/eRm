\encoding{UTF-8}
\name{llraDat2}
\alias{llraDat2}
\docType{data}
\title{An Artificial LLRA Data Set}
\description{Artificial data set of 70 subjects with 5 items, 4 time points and 3 groups for LLRA.}
\usage{llraDat2}
\format{
  A data frame with 70 observations of 21 variables.
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
    \item{\code{group}}{The group membership}
  }
}
\details{
This is a data set as described in Hatzinger and Rusch (2009). 5 items
were measured at 4 time points (in columns). Each persons answers to the
items are recorded in the rows. There are 2
treatment groups and a control group. Treatment group 2 has size, 10,
treatment group 1 has size 20 and the control group has size 40. Item 1 is dichotomous, all others
are polytomous. Item 2, 3, 4 and 5 have 3, 4, 5, 6 categories
respectively.   
}
\references{
Hatzinger, R. and Rusch, T. (2009) IRT models with relaxed assumptions
in eRm: A manual-like instruction. \emph{Psychology Science Quarterly}, \bold{51},
pp. 87--120,
\url{http://erm.r-forge.r-project.org/psq_1_2009_06_87-120.pdf}
}
\examples{
llraDat2
}
\keyword{datasets}
