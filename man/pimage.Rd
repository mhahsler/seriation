\name{pimage}
\alias{pimage}
\alias{pimage.matrix}
\alias{pimage.dist}
\title{Permutation Image Plot}
\description{
  Provides methods for plotting image plots for 
  matrix and \code{dist} objects given a permutation. 
  By default, no permutation is performed.
  This plot can also
  be used as a more versatile replacement of \code{image} plot 
  in \pkg{graphics} based on \pkg{grid}.
}
\usage{
pimage(x, order = NULL, col = NULL, main ="", xlab = "", ylab = "", 
  axes = "auto", zlim=NULL, key=TRUE, key.lab="", symkey=TRUE, 
  upper.tri = TRUE, lower.tri = TRUE, prop = NULL, 
  ..., newpage=TRUE, pop=TRUE, gp=NULL)
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
  \item{key}{ logical; add a color key?}
  \item{key.lab}{ string plotted next to the color key.}
  \item{symkey}{ logical; if \code{x} contains negative values, should the 
    color palate be symmetric (zero is in the middle)>}
  \item{upper.tri, lower.tri}{ a logical indicating whether to show the
      upper or lower triangle of the distance matrix. }
  \item{prop}{logical; draw the cells in the image proportional (defaults to 
    \code{TRUE} for \code{dist} and \code{FALSE} for matrix).}
  \item{\dots}{further arguments passed on to \code{image} in \pkg{graphics}.}
  \item{newpage, pop}{two logical. Start plot on a new page and pop the viewports after plotting (see \pkg{grid}).}
  \item{gp}{a \code{gpar} object (see \pkg{grid}).}
}
\details{
  Plots a matrix in its original row and column orientation. This means, in a
  plot the columns become the x-coordinates and the reversed rows the
  y-coordinates.

  If \code{x} is of class \code{dist} it is converted to full-storage
  representation before plotting.
  
  The viewports used for plotting are called: \code{"plot"}, \code{"image"}
  and \code{"colorkey"}.
}
\seealso{
\code{\link{seriate}}, 
\code{\link{hmap}},
\code{\link{dissplot}} and
\code{\link[graphics]{image}}.}
\author{Christian Buchta and Michael Hahsler}
\examples{
x <- matrix(sample(c(FALSE, TRUE), 300, rep=TRUE), ncol=10, 
  dimnames = list(1:30, LETTERS[1:10]))

## matrix (large values are dark/black)
pimage(x, main = "Random data", key = FALSE)

## plot seriated matrix (use red, proportional display and plot all axes)
pimage(x, seriate(x), col = c("white", "red"), 
  prop = TRUE, axes="both", main = "Reordered data", key = FALSE)

## show correlation (for neg. values a diverging color scheme is 
## used automatically) 
pimage(cor(x), prop=TRUE)
  
## distances (note that low distances are represented dark!)
d <- dist(x, method = "binary")
pimage(d, upper.tri = FALSE, main = "Distances")

pimage(d, seriate(d), upper.tri = FALSE, main = "Distances", axes = "both")

## add to the plot using grid (use pop = FALSE) 
library(grid)
pimage(x, pop = FALSE)
downViewport(name = "image")
  
## highlight cell 7/5 with a red arrow
grid.lines(x = c(5, 7), y = c(3, 5), arrow = arrow(), 
  default.units = "native", gp = gpar(col="red", lwd = 3))
  
## add a red box around rows 15 and 16  
grid.rect(x = 0.5, y = 15.5, width = ncol(x), height = 2, 
  just = "left",
  default.units = "native", gp = gpar(col="red", lwd = 3, fill = NA))

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
}
\keyword{hplot}
