\encoding{UTF-8}
\name{plotPImap}
\alias{plotPImap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Person-Item Map}
\description{
    A person-item map displays the location of item (and threshold) parameters
    as well as the distribution of person parameters.along the latent dimension.
    Person-item maps are useful to compare the range and position of the item measure distribution
    (lower panel) to the range and position of the person measure distribution (upper panel).
    Items should ideally be located along the whole scale to meaningfully measure
    the `ability' of all persons.
}
\usage{
plotPImap(object, item.subset = "all", sorted = FALSE,
   main = "Person-Item Map", latdim = "Latent Dimension",
   pplabel = "Person\nParameter\nDistribution", cex.gen = 0.7,
   xrange = NULL, warn.ord = TRUE, warn.ord.colour = "black",
   irug = TRUE, pp = NULL)
}
\arguments{
  \item{object}{Object of class \code{Rm} or \code{dRm}}
  \item{item.subset}{Subset of items to be plotted. Either a numeric vector indicating
       the column in \code{X} or a character vector indicating the column name.
       If \code{"all"}, all items are plotted. The number of items to be plotted must be > 1.}
  \item{sorted}{ If \code{TRUE}, the items are sorted in increasing order according to their location
       on the latent dimension.}
  \item{main}{Main title of the plot.}
  \item{latdim}{Label of the x-axis, i.e., the latent dimension.}
  \item{pplabel}{Title for the upper panel displaying the person parameter distribution}
  \item{cex.gen}{\code{cex} as a graphical parameter
       specifies a numerical value giving the amount by which plotting text and symbols should be
       magnified relative to the default. Here \code{cex.gen} applies to all text labels. The default is 0.7.}
  \item{xrange}{Range for the x-axis}
  \item{warn.ord}{If \code{TRUE} (the default) asterisks are displayed in the right margin of the lower
       panel to indicate nonordinal threshold locations for polytomous items.}
  \item{warn.ord.colour}{Nonordinal threshold locations for polytomous
    items are coloured with this colour to make them more visible.  This
    is especially useful when there are many items so that the plot is
    quite dense.  The default is \code{"black"}, so that there is no
    distinction made.}
  \item{irug}{If \code{TRUE} (the default), all thresholds are plotted below the person distribution
       to indicate where the included items are most informative.}
  \item{pp}{If non-\code{NULL}, this contains the
    \code{person.parameter} data of the data object, avoiding the
    need to recalculate it.}
}
\details{
  Item locations are displayed with bullets, threshold locations with circles.
}
%\value{}
\references{Bond, T.G., and Fox Ch.M. (2007) Applying the Rasch Model. Fundamental Measurement in the Human Sciences.
2nd Edition. Lawrence Erlbaum Associates.
}
\author{Patrick Mair, Reinhold Hatzinger, patches from Julian Gilbey and Marco J. Maier}
%\note{}
%\seealso{}
\examples{
res <- PCM(pcmdat)
plotPImap(res, sorted=TRUE)
}
\keyword{models}
