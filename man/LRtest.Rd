\encoding{UTF-8}
\name{LRtest}
\alias{LRtest.Rm}
\alias{LRtest}
\alias{print.LR}
\alias{summary.LR}
\alias{plotGOF}
\alias{plotGOF.LR}
\title{Computation of Andersen's LR-test.}
\description{This LR-test is based on subject subgroup splitting.}
\usage{\method{LRtest}{Rm}(object, splitcr = "median", se = TRUE)

\method{plotGOF}{LR}(x, beta.subset = "all", main = "Graphical Model Check", xlab, ylab,
    tlab = "item", xlim, ylim, type = "p", pos = 4, conf = NULL, ctrline = NULL, 
    smooline = NULL, asp = 1, x_axis = TRUE, y_axis = TRUE, set_par = TRUE, 
    reset_par = TRUE, \dots)}
\arguments{
  \item{object}{Object of class \code{"Rm"}.}
  \item{splitcr}{Split criterion for subject raw score splitting.
    \code{"all.r"} corresponds to a full raw score split, \code{"median"} uses the median as split criterion, \code{"mean"} performs a mean split.
    Optionally \code{splitcr} can also be a vector which assigns each person to a certain subgroup (e.g., following an external criterion).
    This vector can be numeric, character or a factor.}
  \item{se}{controls computation of standard errors in the submodels (default: \code{TRUE}).}

  \item{x}{Object of class \code{"LR"}. Also used for visualizing the fit of single items.}
  \item{beta.subset}{If \code{"all"}, all items are plotted. Otherwise numeric subset vector can be specified.}
  \item{main}{Title of the plot.}
  \item{xlab}{Label on \eqn{x}{x}-axis, default gives name of \code{splitcr} and level.}
  \item{ylab}{Label on \eqn{y}{y}-axis, default gives name of \code{splitcr} and level.}
  \item{tlab}{Specification of item labels: \code{"item"} prints the item names, \code{"number"} gives integers corresponding to order of the beta parameters, if \code{"none"} no labels are printed.
    \code{"identify"} allows for an interactive labelling.
    Initially no labels are printed, after clicking close to an item point the corresponding label is added.
    The identification process is terminated by clicking the second button and selecting 'Stop' from the menu, or from the 'Stop' menu on the graphics window.
    For more information and basic operation see \code{\link{identify}}.}
  \item{xlim}{Limits on \eqn{x}{x}-axis.}
  \item{ylim}{Limits on \eqn{y}{y}-axis.}
  \item{type}{Plotting type (see \code{\link[graphics]{plot}}).}
  \item{pos}{Position of the item label (see \code{\link[graphics]{text}}).}
  \item{conf}{for plotting confidence ellipses for the item parameters.
    If \code{conf = NULL} (the default) no ellipses are drawn.
    Otherwise, \code{conf} must be specified as a list with optional elements: \code{gamma}, is the confidence level (numeric), \code{col} and \code{lty}, color and linetype (see \code{\link[graphics:par]{par}}), \code{which} (numeric index vector) specifying for which items ellipses are drawn (must be a subset of \code{beta.subset}), and \code{ia}, logical, if the ellipses are to be drawn interactively (cf., \code{tlab = "identify"} above).
    For details about the default behavior, if \code{conf} is specified as a an empty list, see Details and Examples below.
    To use \code{conf}, the LR object \code{x} has to be generated using the option \code{se = TRUE} in \code{LRtest()}.
    For specification of \code{col} and \code{which} see Details and Examples below.}
  \item{ctrline}{for plotting confidence bands (control lines, cf. eg. Wright and Stone, 1999).
    If \code{ctrline = NULL} (the default) no lines are drawn.
    Otherwise, \code{ctrline} must be specified as a list with optional elements: \code{gamma}, is the confidence level (numeric), \code{col} and \code{lty}, color and linetype (see \code{\link[graphics:par]{par}}).
    If \code{ctrline} is specified as \code{ctrline = list()}, the default values \code{conf = list(gamma = 0.95, col = "blue", lty = "solid")} will be used.
    See examples below.
    To use \code{ctrline}, the LR object \code{x} has to be generated using the option \code{se = TRUE} in \code{LRtest()}.}
  \item{smooline}{spline smoothed confidence bands; must be specified as a list with optional elements: \code{gamma}, is the confidence level (numeric), \code{col} and \code{lty}, color and linetype, \code{spar} as smoothing parameter (see \code{\link{smooth.spline}}).}  
  \item{asp}{sets the \eqn{y/x}{y/x} ratio of the plot (see \code{\link[graphics]{plot.window}}).}
  \item{x_axis}{if \code{TRUE}, the \eqn{x}{x}-axis will be plotted.}
  \item{y_axis}{if \code{TRUE}, the \eqn{y}{y}-axis will be plotted.}
  \item{set_par}{if \code{TRUE}, graphical parameters will be set by the function to optimize the plot's appearance. Unless \code{reset_par = FALSE}, these will be reset to the previous \code{\link[graphics:par]{par}} settings.}
  \item{reset_par}{if \code{TRUE}, graphical parameters will be reset to defaults via \code{\link[graphics:par]{par}()} after plotting (only if \code{set_par = TRUE}). To make adjustments \emph{after} using \code{plotGOF}, this reset can be switched off. Note that the changed graphical parameters will remain in place unless they are redefined (using \code{\link[graphics:par]{par}()}) or the device is closed.}
  \item{...}{additional parameters.}
}
\details{
If the data set contains missing values and \code{mean} or \code{median} is specified as split criterion, means or medians are calculated for each missing value subgroup and consequently used for raw score splitting.

When using interactive selection for both labelling of single points (\code{tlab = "identify"} and drawing confidence ellipses at certain points (\code{ia = TRUE}) then first all plotted points are labelled and afterwards all ellipses are generated.
Both identification processes can be terminated by clicking the second (right) mouse button and selecting `Stop' from the menu, or from the `Stop' menu on the graphics window.

Using the specification \code{which} in allows for selectively drawing ellipses for certain items only, e.g., \code{which = 1:3} draws ellipses for items 1 to 3 (as long as they are included in \code{beta.subset}).
The default is drawing ellipses for all items.
The element \code{col} in the \code{conf} list can either be a single color specification such as \code{"blue"} or a vector with color specifications for all items.
The length must be the same as the number of ellipses to be drawn.
For color specification a palette can be set up using standard palettes (e.g., \code{\link{rainbow}}) or palettes from the \code{colorspace} or \code{RColorBrewer} package.
An example is given below.

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

Mair, P., and Hatzinger, R. (2007). Extended Rasch modeling: The \pkg{eRm} package for the application of IRT models in R. Journal of Statistical Software, 20(9), 1-20.

Mair, P., and Hatzinger, R. (2007). CML based estimation of extended Rasch models with the \pkg{eRm} package in R. Psychology Science, 49, 26-43.

Wright, B.D., and Stone, M.H. (1999). Measurement essentials. Wide Range Inc., Wilmington. (\url{https://www.rasch.org/measess/me-all.pdf} 28Mb).
}
\author{Patrick Mair, Reinhold Hatzinger, Marco J. Maier, Adrian Bruegger}
\seealso{\code{\link{Waldtest}}, \code{\link{MLoef}}}
\examples{
# the object used is the result of running ... RM(raschdat1)
res <- raschdat1_RM_fitted       # see ? raschdat1_RM_fitted

# LR-test on dichotomous Rasch model with user-defined split
splitvec <- sample(1:2, 100, replace = TRUE)
lrres <- LRtest(res, splitcr = splitvec)
lrres
summary(lrres)

\dontrun{
# goodness-of-fit plot with interactive labelling of items w/o standard errors
plotGOF(lrres, tlab = "identify")}

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

\dontrun{
# goodness-of-fit plot for items 1, 14, 24, and 25
# for items 1 and 24 additional 95 percent confidence ellipses
# using colors for these 2 items from the colorspace package
library("colorspace")
my_colors <- rainbow_hcl(2)
plotGOF(lrres2, beta.subset = c(14, 25, 24, 1),
        conf = list(which = c(1, 14), col = my_colors))}

# first, save current graphical parameters in an object
old_par <- par(mfrow = c(1, 2), no.readonly = TRUE)
# plots
plotGOF(lrres2, ctrline = list(gamma = 0.95, col = "red", lty = "dashed"),
  xlim = c(-3, 3), x_axis = FALSE, set_par = FALSE)
axis(1, seq(-3, 3, .5))

plotGOF(lrres2, conf = list(), xlim = c(-3, 3), x_axis = FALSE, set_par = FALSE)
axis(1, seq(-3, 3, .5))
text(-2, 2, labels = "Annotation")
# reset graphical parameters
par(old_par)
}
\keyword{models}
