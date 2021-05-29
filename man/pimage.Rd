\name{pimage}
\alias{pimage}
\alias{pimage.matrix}
\alias{pimage.dist}
\alias{ggpimage}
\alias{ggpimage.matrix}
\alias{ggpimage.dist}
\title{Permutation Image Plot}
\description{
  Provides methods for plotting image plots for
  matrix and \code{dist} objects given a permutation.
  By default, no permutation is performed.
  This plot represents the rows in the order in which they apprear in the matrix
  and provides a color key.
  Implementations based on the \pkg{grid} gaphics engine and based n \pkg{ggplot2} are provided.
}
\usage{
# grid-based image plot
pimage(x, order = NULL, col = NULL, main ="", xlab = "", ylab = "",
  axes = "auto", zlim = NULL, key = TRUE, key.lab = "", symkey = TRUE,
  upper.tri = TRUE, lower.tri = TRUE, interpolate = FALSE, prop = NULL,
  ..., newpage = TRUE, pop = TRUE, gp = NULL)

# ggplot2-based image plot
ggpimage(x, order = NULL, axes = "auto",
  upper.tri = TRUE, lower.tri = TRUE)
}

\arguments{
  \item{x}{a matrix or an object of class \code{dist}.}
  \item{order}{an object of class \code{ser_permutation}. If \code{NULL} the order
      in \code{x} is plotted.}
  \item{col}{a list of colors used. If \code{NULL}, a gray scale is used (for
      matrix larger values are displayed darker and for \code{dist} smaller
      distances are darker). For matrices containing logical data, black and
      white is used. For matrices containing negative values a symmetric
      diverging color palette is used.}
  \item{main}{plot title.}
  \item{xlab, ylab}{labels for the x and y axes.}
  \item{axes}{ a character string indicating if axes labels (column and row
    names of \code{x}) should be potted. Possible values are
    \code{"auto"} (only plot if less then 25 labels),
    \code{"x"},
    \code{"y"},
    \code{"both"} and
    \code{"none"}.}
  \item{zlim}{ vector with two elements giving the range (min, max) for
      representing the values in the matrix. }
  \item{key}{ logical; add a color key? No key is available for logical
    matrices. }
  \item{key.lab}{ string plotted next to the color key.}
  \item{symkey}{ logical; if \code{x} contains negative values, should the
    color palate be symmetric (zero is in the middle)>}
  \item{upper.tri, lower.tri}{ a logical indicating whether to show the
      upper or lower triangle of the distance matrix. }
  \item{interpolate}{logical; interpolate color between matrix cells. }
  \item{prop}{logical; draw the cells in the image proportional (defaults to
    \code{TRUE} for \code{dist} and \code{FALSE} for matrix).}
  \item{\dots}{further arguments passed on to \code{image} in \pkg{graphics}.}
  \item{newpage, pop}{two logical. Start plot on a new page and pop the viewports after plotting (see \pkg{grid}).}
  \item{gp}{a \code{gpar} object (see \pkg{grid}).}
}
\details{
  Plots a matrix in its original row and column orientation. This means, in a
  plot the columns become the x-coordinates and the rows the
  y-coordinates (in reverse order).

  If \code{x} is of class \code{dist} it is converted to full-storage
  representation before plotting.

  \bold{Grid-based plot:}
  The viewports used for plotting are called: \code{"plot"}, \code{"image"}
  and \code{"colorkey"}.
  \emph{Note:} Since \code{pimage} uses \pkg{grid}, it should not be mixed with base R primitive plotting
  functions, but the appropriate functions in \code{\link[grid]{grid-package}}.

  \bold{ggplot2-based plot:} A ggplot2 object is returned. Colors, axis limites and other visual aspects
  can be added using standard ggplot2 functions
  (\code{labs}, \code{scale_fill_continuous}, \code{labs}, etc.).
}
\seealso{
\code{\link[grid]{grid-package}},
\code{\link{seriate}},
\code{\link{hmap}},
\code{\link{dissplot}} and
\code{\link[graphics]{image}}.}
\author{Christian Buchta and Michael Hahsler}
\examples{
x <- matrix(sample(c(FALSE, TRUE), 300, rep = TRUE), ncol = 10,
  dimnames = list(1:30, LETTERS[1:10]))

## matrix (TRUE values are black)
pimage(x, main = "Random data", key = FALSE)

## plot seriated matrix (use red, proportional display and plot all axes)
pimage(x, seriate(x), col = c("white", "red"),
  prop = TRUE, axes = "both", main = "Reordered data", key = FALSE)

## show correlation (with neg. values, a diverging color scheme is
## used automatically)
corr <- cor(x)
pimage(corr, order = seriate(corr), prop = TRUE)

## distances (note that low distances are represented dark!)
d <- dist(x, method = "binary")
pimage(d, upper.tri = FALSE, main = "Distances")

pimage(d, order = seriate(d), upper.tri = FALSE, main = "Distances", axes = "both")

## add to the plot using functions in package grid
library("grid")
pimage(x, pop = FALSE)           ### pop = FALSE allows us to manipulate viewports
downViewport(name = "image")

## highlight cell 7/5 with a red arrow
grid.lines(x = c(5, 7), y = c(3, 5), arrow = arrow(),
  default.units = "native", gp = gpar(col = "red", lwd = 3))

## add a red box around rows 15 and 16
grid.rect(x = 0.5, y = 15.5, width = ncol(x), height = 2,
  just = "left",
  default.units = "native", gp = gpar(col = "red", lwd = 3, fill = NA))

## remove the viewports
popViewport(0)

## put several pimages on a page (uses viewports and newpage = FALSE)
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow = 1, ncol = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

## seriate matrix
o <- seriate(x)
pimage(x, o, main = "Random data", prop = TRUE, axes = "both", key = FALSE,
  newpage = FALSE)

upViewport(1)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))

## add the reordered disimilarity matrix for rowa
pimage(d, o[[1]], upper.tri = FALSE, main = "Distances", axes = "both",
  key = FALSE, newpage = FALSE)

upViewport(1)
popViewport(0)

## ggplot2 Examples
if (require("ggplot2")) {
  x <- matrix(sample(c(FALSE, TRUE), 300, rep = TRUE), ncol = 10,
    dimnames = list(1:30, LETTERS[1:10]))

  ## matrix (TRUE values are dark). We also hide the legend
  ggpimage(x, axes = "both") +
    labs(title = "Random data") +
    theme(legend.position = "none")

  ## show reordered correlation (with neg. values, a diverging color scheme is
  ## used automatically)
  corr <- cor(x)
  ggpimage(corr, order = seriate(corr)) +
    theme(aspect.ratio = 1)

  ## distances (distances are shown with an aspect ratio of 1 and dark means low distance)
  ## the suppressed triangle is set to white with na.value = "white"
  d <- dist(x, method = "binary")
  ggpimage(d, order = seriate(d), upper.tri = FALSE, axes = "both") + labs(title = "Distances")

  ## replace default colors using Colorbrewer and colorspace
  ggpimage(d, order = seriate(d), upper.tri = FALSE) +
    scale_fill_distiller(palette = "Spectral", direction = +1)

  library(colorspace)
  ggpimage(d, order = seriate(d), upper.tri = FALSE) +
    scale_fill_continuous_sequential("Reds", rev = FALSE)
}
}
\keyword{hplot}
