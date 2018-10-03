\encoding{UTF-8}
\name{plotDIF}
\alias{plotDIF}
\title{
Confidence intervals plot of item parameter estimates.
}
\description{
Performs an plot of item parameter conficence intervals based on \code{LRtest} subgroup splitting.
}
\usage{
plotDIF(object, item.subset = NULL, gamma = 0.95, main = NULL,
          xlim = NULL,  xlab = " ", ylab=" ", col = NULL,
          distance, splitnames=NULL, leg = FALSE, legpos="bottomleft", ...)
}

\arguments{
  \item{object}{
An object of class \code{LR} (if more objects should be plotted, the argument has to be defined as a \code{list}).
}
  \item{item.subset}{
Subset of items to be plotted. Either a numeric vector indicating the items or a character vector indicating the itemnames.
If nothing is defined (default), all items are plotted.
}
  \item{gamma}{
The level for the item parameter's confidence limits (default is gamma = 0.95).
}
  \item{main}{
Main title for the plot.
}
  \item{xlim}{
Numeric vector of length 2, giving the x coordinates ranges of the plot (the y coordinates depend on the number of
depicted items).
}
  \item{xlab}{
Label for the x axis.
}
  \item{ylab}{
Label for the y axis.
}
  \item{col}{
By default the color for the drawn confidence lines is determined automatically whereas every group (split criterion)
is depicted in the same color.
}
  \item{distance}{
Distance between each item's confidence lines -- if omitted, the distance shrinks with increasing numbers of split criteria. Can be overriden using values in (0, 0.5).
%Distance between the drawn confidence lines, default is division by factor 10 .
}
  \item{splitnames}{
For labeling the splitobjects in the legend (returns a nicer output).
}
  \item{leg}{
If \code{TRUE} a legend is provided by default.
}
  \item{legpos}{
Position of the legend with possible values  \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"},
\code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"} and \code{"center"}. The default value for the legend
is \code{"bottomleft"}.
}
  \item{...}{
Further options to be passed to \code{plot}.
}

}

\details{
 If there are items that cannot be estimated for some reasons, certainly these ones are not plotted.
 For plotting several objects of class \code{LR}, the subgroup splitting by \code{LRtest} has to
 be carried out for the same data set (or at least item subsets of it).

 Plotting a certain subset of items could be useful if the objects of class \code{LR} contain a huge number
 of estimated items.

 The default level for the conficence limits is gamma = 0.95. (If the conficence limits should be
 corrected it is useful to use a  correction, e.g., Bonferroni: 1 - (1 - gamma) / number of estimated items.)
}
\value{
\code{plotCI} returns a list containing the confidence limits of each group in each \code{LRtest} object.
}

\author{Kathrin Gruber, Reinhold Hatzinger}
\seealso{\code{\link{LRtest}}, \code{\link{confint.threshold}}, \code{\link{thresholds}}}
\examples{
# the object used is the result of running   RM(raschdat1)
res <- raschdat1_RM_fitted     # see ? raschdat1_RM_fitted

\dontrun{
# LR-test on dichotomous Rasch model with user-defined split
splitvec <- rep(1:2, each = 50)
lrres <- LRtest(res, splitcr = splitvec)

# LR-test with mean split
lrres2 <- LRtest(res, split = "mean")

# combination of LRtest-objects in a list
RMplotCI <- list(lrres, lrres2)}

# the object raschdat1_RM_plotDIF is the result of the computations outlined
# above and is loaded to save computation time. see ?raschdat1_RM_plotDIF
RMplotCI <- raschdat1_RM_plotDIF

# Confidence intervals plot with default assumptions
plotDIF(RMplotCI)

# Confidence intervals plot with Bonferroni correction
plotDIF(RMplotCI, gamma = (1 - (0.05/10)))

# Confidence intervals plot for an item subset
plotDIF(RMplotCI, item.subset = 1:6)

# with user defined group color and legend
plotDIF(RMplotCI, col = c("red", "blue"), leg = TRUE, legpos = "bottomright")

# with names for the splitobjects
plotDIF(RMplotCI, col = c("red", "blue"), leg = TRUE, legpos = "bottomright",
        splitnames = c(paste("User", 1:2), paste(rep("Mean", 2), 1:2)))
}

\keyword{models}
