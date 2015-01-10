\encoding{UTF-8}
\name{plotPWmap}
\alias{plotPWmap}
\title{Pathway Map}
\description{
    A Bond-and-Fox Pathway Map displays the location of each item or
    each person against its infit t-statistic.  Pathway maps are useful
    for identifying misfitting items or misfitting persons.  Items or
    people should ideally have a infit t-statistic lying between about
    -2 and +2, and these values are marked.
}
\usage{
plotPWmap(object, pmap = FALSE, imap=TRUE,
                 item.subset = "all", person.subset = "all",
                 mainitem = "Item Map", mainperson = "Person Map",
                 mainboth="Item/Person Map",
                 latdim = "Latent Dimension",
                 tlab = "Infit t statistic",
                 pp = NULL, cex.gen = 0.6, cex.pch=1,
                 person.pch = 1, item.pch = 16,
                 personCI = NULL, itemCI = NULL, horiz=FALSE)
}
\arguments{
  \item{object}{Object of class \code{Rm} or \code{dRm}}
  \item{pmap}{Plot a person map if \code{TRUE}; the default is
    \code{FALSE}.}
  \item{imap}{Plot an item map if \code{TRUE} (the default); do not plot
    if \code{FALSE}.  At least one of \code{pmap} and \code{imap} must
    be \code{TRUE}.}
  \item{item.subset}{Subset of items to be plotted for an item map.
    Either a numeric vector indicating the item numbers or a character
    vector indicating the item names.  If \code{"all"}, all items are
    plotted. The number of items to be plotted must be > 1.}
  \item{person.subset}{Subset of persons to be plotted for a person map.
    Either a numeric vector indicating the person numbers or a character
    vector indicating the person names.  If \code{"all"}, all persons are
    plotted. The number of persons to be plotted must be > 1.}
  \item{mainitem}{Main title of an item plot.}
  \item{mainperson}{Main title of a person plot.}
  \item{mainboth}{Main title of a person/item joint plot.}
  \item{latdim}{Label of the y-axis, i.e., the latent dimension.}
  \item{tlab}{Label of the x-axis, i.e., the t-statistic dimension.}
  \item{pp}{If non-\code{NULL}, this contains the
    \code{person.parameter} data of the data object, avoiding the
    need to recalculate it.}
  \item{cex.gen}{\code{cex} as a graphical parameter
    specifies a numerical value giving the amount by which plotting
    text and symbols should be magnified relative to the
    default. Here \code{cex.gen} applies to all text labels. The
    default is 0.6.}
  \item{cex.pch}{applies to all plotting symbols. The
    default is 1.}
  \item{person.pch, item.pch}{Specifies the symbol used for plotting
    person data and item data respectively; the defaults are 1 and 16
    respectively.  See \code{\link{points}} for more information
    about \code{pch} values.}
  \item{personCI, itemCI}{Plotting confidence intervals for the the
    person abilities and item difficulties.  If \code{personCI=NULL}
    (the default) no confidence intervals are drawn for person
    abilities.  Otherwise, specifying \code{personCI} draws
    approximate confidence intervals for each person's ability.
    \code{personCI} must be specified as a list, and the optional
    elements of this list are \code{gamma}, the confidence level,
    \code{col}, colour, and \code{lty}, line type.  If \code{personCI}
    is specified as an empty list, or not all of the list items are
    specified, the default values
    \code{personCI=list(gamma=0.95,col="orange",lty="dotted")} will be
    used.

    The same goes for \code{itemCI}, except that the default settings
    are \code{itemCI=list(gamma=0.95,col="red",lty="dotted")}.}
  \item{horiz}{if \code{TRUE}, the plot is horizontal, i.e., the latent
    dimension is on the x-axis. The default is \code{FALSE}.}
}
\details{
  This code uses vertical(horizontal) error bars rather than circles or boxes to
  indicate standard errors.  It also offers the possibility of plotting
  item or person data on its own; this can considerably simplify the
  reading of the plots for large datasets.
}
%\value{}
\references{
  Bond T.G., Fox C.M. (2007) \emph{Applying the Rasch Model: Fundamental Measurement in the Human Sciences} (2nd ed.) chapter 3, Lawrence Erlbaum Associates, Inc.

  Linacre J.M., Wright B.D. (1994) Dichotomous Infit and Outfit Mean-Square Fit Statistics / Chi-Square Fit Statistics. \emph{Rasch Measurement Transactions} \bold{8:2} p. 350, \url{http://www.rasch.org/rmt/rmt82a.htm}

  Linacre J.M. (2002) What do Infit and Outfit, Mean-square and Standardized mean? \emph{Rasch Measurement Transactions} \bold{16:2} p. 878, \url{http://www.rasch.org/rmt/rmt162f.htm}

  Wright B.D., Masters G.N. (1990) Computation of OUTFIT and INFIT Statistics. \emph{Rasch Measurement Transactions} \bold{3:4} p. 84--85, \url{http://www.rasch.org/rmt/rmt34e.htm}
}
\author{Julian Gilbey}
%\note{}
%\seealso{}
\examples{
res <- PCM(pcmdat)
pparm <- person.parameter(res)
plotPWmap(res, pp = pparm)
plotPWmap(res, pp = pparm, pmap = TRUE)
}
\keyword{models}
