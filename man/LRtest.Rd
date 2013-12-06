\name{LRtest}
\alias{LRtest.Rm}
\alias{LRtest}
\alias{print.LR}
\alias{summary.LR}
\alias{plotGOF}
\alias{plotGOF.LR}
\title{Computation of Andersen's LR-test.}
\description{This LR-test is based on subject subgroup splitting.
}
\usage{
\method{LRtest}{Rm}(object, splitcr = "median", se = TRUE)
\method{plotGOF}{LR}(x, beta.subset = "all", main="Graphical Model Check",
   xlab = NULL, ylab = NULL, tlab = "item",
   ylim = c(-3, 3), xlim = c(-3, 3), type = "p", pos = "4",
   conf = NULL, ctrline = NULL, ...)
%\method{print}{LR}(x,...)
%\method{summary}{LR}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{Rm}.}
  \item{splitcr}{Split criterion for subject raw score splitting. \code{all.r} corresponds to a
  full raw score split, \code{median} uses the median as split criterion, \code{mean} performs a mean-split.
  Optionally \code{splitcr} can also be a vector which assigns each person to a
  certain subgroup (e.g., following an external criterion). This vector can be numeric, character or a factor.}
  \item{se}{If \code{TRUE} standard errors for beta's are computed.}

%Arguments for \code{plotGOF}:
  \item{x}{Object of class \code{LR}. Also used for visualizing the fit of single items.}
  \item{beta.subset}{If \code{"all"}, all items are plotted. Otherwise numeric subset vector can be specified.}
  \item{main}{Main title of the plot.}
  \item{xlab}{Label on x-axis, default gives name of \code{splitcr} and level.}
  \item{ylab}{Label on y-axis, default gives name of \code{splitcr} and level.}
  \item{tlab}{Specification of item labels: \code{"item"} prints the item names, \code{"number"} gives integers
       corresponding to order of the beta parameters, if \code{"none"} no labels are printed.
       \code{"identify"} allows for an interactive labelling. Initially no labels are printed, after clicking
       close to an item point the corresponding label is added. The identification process is terminated by clicking
       the second button and selecting 'Stop' from the menu, or from the 'Stop' menu on the graphics window.
       For more information and basic operation see \code{\link{identify}}.
       }
  \item{xlim}{Limits on x-axis.}
  \item{ylim}{Limits on y-axis.}
  \item{type}{Plotting type.(see \code{\link{plot}})}
  \item{pos}{Position of the item label (see \code{\link{text}})}
  \item{conf}{for plotting confidence ellipses for the item parameters. If \code{conf=NULL}
             (the default) no ellipses are drawn. Otherwise, \code{conf} must be
             specified as a list with optional elements: \code{gamma}, is
             the confidence level (numeric), \code{col} and \code{lty},
             colour and linetype (see \code{\link{par}}), \code{which} (numeric index vector) specifying for which
             items ellipses are drawn (must be a subset of \code{beta.subset}),
             and \code{ia}, logical, if the ellipses are to be drawn interactively (cf.
             \code{tlab="identify"} above). If \code{conf} is specified as a an empty list, %\code{conf=list()},
             the default values \code{conf=list(gamma=0.95, col="red", lty="dashed", ia=FALSE)}
             will be used. See example below. To use \code{conf}, the LR object \code{x} has
             to be generated using the option \code{se=TRUE} in \code{LRtest()}.
             For specification of \code{col} and \code{which} see Details and Examples below.}
  \item{ctrline}{for plotting confidence bands (control lines, cf.eg.Wright and Stone, 1999).
             If \code{ctrline=NULL}
             (the default) no lines are drawn. Otherwise, \code{ctrline} must be
             specified as a list with optional elements: \code{gamma}, is
             the confidence level (numeric), \code{col} and \code{lty}, colour and linetype (see \code{\link{par}}).
             If \code{ctrline} is specified as \code{ctrline=list()},
             the default values \code{conf=list(gamma=0.95, col="blue", lty="solid")}
             will be used. See examples below. To use \code{ctrline}, the LR object \code{x} has
             to be generated using the option \code{se=TRUE} in \code{LRtest()}.
             }
  \item{...}{Additional parameters.}

}
\details{If the data set contains missing values and \code{mean} or \code{median} is specified as splitcriterion,
         means or medians are calculated for each missing value subgroup and consequently used for raw score splitting.

         When using interactive selection for both labelling of single points (\code{tlab = "identify"} and
         drawing confidence ellipses at certain points (\code{ia = TRUE}) then first all plotted points are labelled
         and afterwards all ellipses are generated. Both identification processes can be terminated
         by clicking the second (right) mouse button and selecting `Stop' from the menu, or from the `Stop'
         menu on the graphics window.

         Using the specification \code{which} in allows for selectively drawing ellipses for
         certain items only, e.g., \code{which=1:3} draws ellipses for items 1 to 3 (as long as they are included
         in \code{beta.subset}). The default is drawing ellipses for all items. 
         The element \code{col} in the \code{conf} list can either be a single colour
         specification such as \code{"blue"} or a vector with colour specifications for all items.
         The length must be the same as the number of ellipses to be drawn. For colour specification
         a palette can be set up using standard palettes (e.g. \code{\link{rainbow}}) or palettes from
         the \code{colorspace} or \code{RColorBrewer} package. An example is given below.
         

         \code{summary} and \code{print} methods are available for objects of class \code{LR}.
}
\value{
\code{LRtest} returns an object of class \code{LR} containing:
  \item{LR}{LR-value.}
  \item{df}{Degrees of freedom of the test statistic.}
  \item{Chisq}{Chi-square value with corresponding df.}
  \item{pvalue}{P-value of the test.}
  \item{likgroup}{Log-likelihood values for the subgroups}
  \item{betalist}{List of beta parameters for the subgroups.}
  \item{selist}{List of standard errors of beta's.}
  \item{etalist}{List of eta parameters for the subgroups.}
  \item{spl.gr}{Names and levels for \code{splitcr}.}
  \item{call}{The matched call.}
  \item{fitobj}{List containing model objects from subgroup fit.}
}
\references{
Fischer, G. H., and Molenaar, I. (1995). Rasch Models - Foundations, Recent Developements, and Applications. Springer.

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The eRm package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the eRm package in R. Psychology Science, 49, 26-43.

Wright, B.D.,  and Stone, M.H. (1999). Measurement essentials. Wide Range Inc., Wilmington.
        (\url{http://www.rasch.org/measess/me-all.pdf} 28Mb).
}

\author{Patrick Mair, Reinhold Hatzinger}
%\note{}
\seealso{\code{\link{Waldtest}}, \code{\link{MLoef}}}
\examples{
# the object used is the result of running   RM(raschdat1)
res <- raschdat1_RM_fitted     # see ? raschdat1_RM_fitted

# LR-test on dichotomous Rasch model with user-defined split
splitvec <- sample(1:2, 100, replace = TRUE)
lrres <- LRtest(res, splitcr = splitvec)
lrres
summary(lrres)


\dontrun{
# goodness-of-fit plot with interactive labelling of items w/o standard errors
plotGOF(lrres, tlab = "identify", se = FALSE)}


# LR-test with a full raw-score split
X <- sim.rasch(1000, -2:2, seed = 5)
res2 <- RM(X)
full_lrt <- LRtest(res2, splitcr = "all.r")
full_lrt

\dontrun{
# LR-test with mean split, standard errors for beta's
lrres2 <- LRtest(res, split = "mean")}
# to save computation time, the results are loaded from raschdat1_RM_lrres2
lrres2 <- raschdat1_RM_lrres2                    # see ?raschdat1_RM_lrres2

# goodness-of-fit plot
# additional 95 percent control line with user specified style
plotGOF(lrres2, ctrline = list(gamma = 0.95, col = "red", lty = "dashed"))

# goodness-of-fit plot for items 1, 14, 24, and 25
# additional 95 percent confidence ellipses, default style
plotGOF(lrres2, beta.subset = c(14, 25, 24, 1), conf = list())

# goodness-of-fit plot for items 1, 14, 24, and 25
# for items 1 and 24 additional 95 percent confidence ellipses
# using colours for these 2 items from the colorspace package
\dontrun{
library(colorspace)
colors <- rainbow_hcl(2)
plotGOF(lrres2, beta.subset = c(14, 25, 24, 1),
        conf = list(which = c(1, 14), col = colors))}
}

\keyword{models}
