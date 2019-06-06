\encoding{UTF-8}
\name{plotICC}
\alias{plotICC}
\alias{plotICC.Rm}
\alias{plotjointICC}
\alias{plotjointICC.dRm}

\title{ICC Plots}
\description{Plot functions for visualizing the item characteristic curves}
\usage{
\method{plotICC}{Rm}(object, item.subset = "all", empICC = NULL, empCI = NULL,
   mplot = NULL, xlim = c(-4, 4), ylim = c(0, 1),
   xlab = "Latent Dimension", ylab = "Probability to Solve", main=NULL,
   col = NULL, lty = 1, legpos = "left", ask = TRUE, ...)
\method{plotjointICC}{dRm}(object, item.subset = "all", legend = TRUE,
   xlim = c(-4, 4), ylim = c(0, 1), xlab = "Latent Dimension",
   ylab = "Probability to Solve", lty = 1, legpos = "topleft",
   main="ICC plot",col=NULL,...)
}
\arguments{
  \item{object}{object of class \code{Rm} or \code{dRm}}
  \item{item.subset}{Subset of items to be plotted. Either a numeric vector indicating
        the column in \code{X} or a character vector indiciating the column name.
        If \code{"all"} (default), all items are plotted.}
  \item{empICC}{Plotting the empirical ICCs for objects of class \code{dRm}.
        If \code{empICC=NULL}
        (the default) the empirical ICC is not drawn. Otherwise, \code{empICC} must be
        specified as a list where the first element must be one of
        \code{"raw"}, \code{"loess"}, \code{"tukey"}, \code{"kernel"}. The other optional elements are
        \code{smooth} (numeric), \code{type} (line type for empirical ICCs,
        useful values are \code{"p"} (default), \code{"l"}, and \code{"b"},
          see graphics parameter \code{type} in \code{\link{plot.default}}),
        \code{pch}, \code{col}, and \code{lty}, plotting `character', colour and linetype
        (see \code{\link{par}}). See details and examples below.
  }
  \item{empCI}{Plotting confidence intervals for the the empirical ICCs.
        If \code{empCI=NULL} (the default) no confidence intervals are drawn.
        Otherwise, by specifying \code{empCI} as a list gives `exact' confidence
        intervals for each point of the empirical ICC.
        The optional elements of this list are \code{gamma}, the confidence level,
        \code{col}, colour, and \code{lty}, line type. If \code{empCI} is specified
        as an empty list,
        the default values \code{empCI=list(gamma=0.95,col="red",lty="dotted")}
        will be used.
  }
  \item{mplot}{if \code{NULL} the default setting is in effect. For models of class \code{dRm} this
        is \code{mplot = TRUE}, i.e.,
        the ICCs for up to 4 items are plotted in one figure. For \code{Rm}
        models the default is \code{FALSE} (each item in one figure) but may be set to \code{TRUE}.
  }
  \item{xlab}{Label of the x-axis.}
  \item{ylab}{Label of the y-axis.}
  \item{xlim}{Range of person parameters.}
  \item{ylim}{Range for probability to solve.}
  \item{legend}{If \code{TRUE}, legend is provided, otherwise the ICCs are labeled.}
  \item{col}{If not specified or \code{NULL}, line colors are determined automatically.
             Otherwise, a scalar or vector with appropriate color specifications may be supplied
             (see \code{\link{par}}).}
  \item{lty}{Line type.}
  \item{main}{Title of the plot.}
  \item{legpos}{Position of the legend with possible values  \code{"bottomright"},
               \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"},
               \code{"topright"}, \code{"right"} and \code{"center"}.
               If \code{FALSE} no legend is displayed.}
  \item{ask}{If \code{TRUE} (the default) and the \code{R} session is interactive the user is asked for input,
             before a new figure is drawn. \code{FALSE} is only useful if automated figure export is
             in effect, e.g., when using \code{\link{Sweave}}.}
  \item{\ldots}{Additional plot parameters.}
}
\details{Empirical ICCs for objects of class \code{dRm} can be plotted using the option \code{empICC}, a
         list where the first element specifies the type of calculation of the empirical values.
         If \code{empICC=list("raw", other specifications)}
         relative frequencies of the positive responses are calculated for each rawscore group and plotted
         at the position of the corresponding person parameter. The other options use the default versions
         of various smoothers: \code{"tukey"} (see \code{\link{smooth}}), \code{"loess"} (see \code{\link{loess}}),
         and \code{"kernel"} (see \code{\link{ksmooth}}). For \code{"loess"} and \code{"kernel"} a further
         element, \code{smooth},
         may be specified to control the span (default is 0.75) or the bandwith (default is 0.5),
         respectively. For example, the specification could be \code{empirical = list("loess", smooth=0.9)}
         or \code{empirical = list("kernel",smooth=2)}.
         Higher values result in smoother estimates of the empirical ICCs.

         The optional confidence intervals are obtained by a procedure first given in
         Clopper and Pearson (1934) based on the beta distribution (see \code{\link{binom.test}}).


}
\note{For most of the plot options see  \code{\link{plot}} and \code{\link{par}}.}
%\value{}
%\references{}
\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{plotGOF}}}
\examples{
\dontrun{
# Rating scale model, ICC plot for all items
rsm.res <- RSM(rsmdat)
thresholds(rsm.res)
plotICC(rsm.res)

# now items 1 to 4 in one figure without legends
plotICC(rsm.res, item.subset = 1:4, mplot = TRUE, legpos = FALSE)

# Rasch model for items 1 to 8 from raschdat1
# empirical ICCs displaying relative frequencies (default settings)
rm8.res <- RM(raschdat1[,1:8])
plotICC(rm8.res, empICC=list("raw"))

# the same but using different plotting styles
plotICC(rm8.res, empICC=list("raw",type="b",col="blue",lty="dotted"))

# kernel-smoothed empirical ICCs using bandwidth = 2
plotICC(rm8.res, empICC = list("kernel",smooth=3))

# raw empirical ICCs with confidence intervals
# displaying only items 2,3,7,8
plotICC(rm8.res, item.subset=c(2,3,7,8), empICC=list("raw"), empCI=list())

# Joint ICC plot for items 2, 6, 8, and 15 for a Rasch model
res <- RM(raschdat1)
plotjointICC(res, item.subset = c(2,6,8,15), legpos = "left")
}
}
\keyword{models}
