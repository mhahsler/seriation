\name{pimage}
\alias{pimage}
\alias{pimage.matrix}
\alias{pimage.dist}
\alias{ggpimage}
\alias{ggpimage.matrix}
\alias{ggpimage.dist}
\title{Permutation Image Plot}
\description{
  Provides methods for matrix shading, i.e., displaying a color image for
  matrix and \code{dist} objects given an optional permutation.
  The plot arranges colored rectangles to
  represent the matrix value. Columns and rows appear in the order in the matrix.
  This visualization is also know asi a heatmap.
  Implementations based on the \pkg{grid} graphics engine and based n \pkg{ggplot2} are provided.
}
\usage{
# grid-based image plot
pimage(x, order = NULL, col = NULL, main ="", xlab = "", ylab = "",
  zlim = NULL, key = TRUE, keylab = "", symkey = TRUE,
  upper_tri = TRUE, lower_tri = TRUE, row_labels = NULL, col_labels = NULL,
  prop = TRUE, flip_axes = FALSE, reverse_columns = FALSE,
  ..., newpage = TRUE, pop = TRUE, gp = NULL)

# ggplot2-based image plot
ggpimage(x, order = NULL,
  upper_tri = TRUE, lower_tri = TRUE, row_labels = NULL, col_labels = NULL,
  prop = TRUE, flip_axes = FALSE, reverse_columns = FALSE)
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
  \item{zlim}{ vector with two elements giving the range (min, max) for
      representing the values in the matrix. }
  \item{key}{ logical; add a color key? No key is available for logical
    matrices. }
  \item{keylab}{ string plotted next to the color key.}
  \item{symkey}{ logical; if \code{x} contains negative values, should the
    color palate be symmetric (zero is in the middle)>}
  \item{upper_tri, lower_tri}{ a logical indicating whether to show the
      upper or lower triangle of the distance matrix. }
  \item{row_labels, col_labels}{ a logical indicating if row and column labels in \code{x} should be displayed.
    If \code{NULL} then labels are displayed if the \code{x} contains the appropriate dimname and the number of labels is 25 or less. A character vector of the appropriate length with labels can also be supplied. }
  \item{prop}{logical; change the aspect ratio so cells in the image have a equal width and height.}
  \item{flip_axes}{logical; exchange rows and columns for plotting.}
  \item{reverse_columns}{logical; revers the order of how the columns are displayed.}
  \item{\dots}{further arguments are ignored. }
  \item{newpage, pop, gp}{Start plot on a new page, pop the viewports after plotting, and
    use the supplied \code{gpar} object (see \pkg{grid}).}
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

  \bold{ggplot2-based plot:} A ggplot2 object is returned. Colors, axis limits and other visual aspects
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

# Matrix for logical values. TRUE values are dark and no color key is shown. There are too many
# Row labels (>25) so they are suppressed.
pimage(x)

# Show all labels and flip axes or reverse columns
pimage(x, row_labels = TRUE, col_labels = TRUE, flip_axes = TRUE)
pimage(x, row_labels = TRUE, col_labels = TRUE, reverse_columns = TRUE)

# Reorder matrix, use custom colors, and add a title.
pimage(x, order = seriate(x), row_labels = TRUE, col_labels = TRUE,
  col = c("white", "red"), main = "Random Data (Reordered)")

# Show reordered correlation (between column vectors). For matrices with negative values,
# a diverging color scheme is  used automatically. The aspect ratio is fixed to 1:1 and only
# the lower triangle is shown.
corr <- cor(x)
pimage(corr, order = seriate(corr), upper_tri = FALSE, main = "Random Data (Correlation)")

# Show a reordered distance matrix (distances between rows).
# Dark means low distance. The aspect ratio is automatically fixed to 1:1.
# The upper triangle is suppressed triangle
d <- dist(x, method = "binary")
pimage(d,  order = seriate(d),
  main = "Random Data (Distances)")

# Show only distances that are smaller than .5 using limits on z.
pimage(d,  order = seriate(d),
  main = "Random Data (Distances + Theshold)", zlim = c(0, .5))

# Add to the plot using functions in package grid
library("grid")
pimage(x, pop = FALSE)           ### pop = FALSE allows us to manipulate viewports
downViewport(name = "image")

# Highlight cell column 7 (G) / row 5 (from top)/col with a red arrow starting at 5/2
# Note: columns are x and rows are y.
grid.lines(x = c(5, 7), y = c(2, 5), arrow = arrow(),
  default.units = "native", gp = gpar(col = "red", lwd = 3))

# add a red box around rows 15 and 16
grid.rect(x = 0.5, y = 15.5, width = ncol(x), height = 2,
  just = "left",
  default.units = "native", gp = gpar(col = "red", lwd = 3, fill = NA))

## remove the viewports
popViewport(0)

## put several pimages on a page (uses viewports and newpage = FALSE)
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2,
  widths = unit(c(.4, .6), unit = "npc"))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

## seriate matrix
o <- seriate(x)
pimage(x, o, labCol = FALSE, main = "Random Data",
  newpage = FALSE)

upViewport(1)
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))

## add the reordered dissimilarity matrix for rows
d <- dist(x)
pimage(d, o[[1]], main = "Distances",
  newpage = FALSE)

upViewport(1)
popViewport(0)

## ggplot2 Examples
if (require("ggplot2")) {
  x <- matrix(sample(c(FALSE, TRUE), 300, rep = TRUE), ncol = 10,
    dimnames = list(1:30, LETTERS[1:10]))

  # Matrix for logical values. TRUE values are dark. There are too many
  # Row labels (>25) so they are suppressed.
  ggpimage(x)

  # Show all labels and flip axes or reverse columns
  ggpimage(x, flip_axes = TRUE, row_labels = TRUE, col_labels = TRUE)
  ggpimage(x, reverse_columns = TRUE, row_labels = TRUE, col_labels = TRUE)

  # Add lines
  ggpimage(x) +
    geom_hline(yintercept = seq(0, nrow(x)) + .5) +
    geom_vline(xintercept = seq(0, ncol(x)) + .5)

  # Reorder matrix, use custom colors, add a title,
  # and hide colorkey.
  ggpimage(x, order = seriate(x), row_labels = TRUE, col_labels = TRUE) +
    scale_fill_manual(values = c("white", "red")) +
    theme(legend.position = "none") +
    labs(title = "Random Data")

  # Show reordered correlation (between column vectors). For matrices with negative values,
  # a diverging color scheme is  used automatically. We show onlt the lower triangle.
  corr <- cor(x)
  ggpimage(corr, order = seriate(corr), upper_tri = FALSE) +
    labs(title = "Random Data", subtitle = "Correlation")

  # Show a reordered distance matrix (distances between rows).
  # Dark means low distance.
  d <- dist(x, method = "binary")
  ggpimage(d, order = seriate(d)) +
    labs(title = "Random Data", subtitle = "Distances")

  # Show also upper triangle
  ggpimage(d, order = seriate(d), upper_tri = TRUE) +
    labs(title = "Random Data", subtitle = "Distances")


  # Show only distances that are smaller than .5 using limits on fill.
  ggpimage(d,  order = seriate(d)) +
    labs(title = "Random Data (Distances + Theshold)") +
    scale_fill_gradient(low = "darkblue", high = "lightgray", limit = c(0, .5), na.value = "white")

  # Use ggplot2 themes with theme_set
  old_theme <- theme_set(theme_linedraw())
  ggpimage(d, order = seriate(d)) +
    labs(title = "Random Data (Distances)")
  theme_set(old_theme)

  # Use custom color palettes: Gray scale, Colorbrewer (provided in ggplot2) and colorspace
  ggpimage(d, order = seriate(d), upper_tri = FALSE) +
    scale_fill_gradient(low = "black", high = "white", na.value = "white")

  ggpimage(d, order = seriate(d), upper_tri = FALSE) +
    scale_fill_distiller(palette = "Spectral", direction = +1, na.value = "white")

  ggpimage(d, order = seriate(d), upper_tri = FALSE) +
    colorspace::scale_fill_continuous_sequential("Reds", rev = FALSE, na.value = "white")
}
}
\keyword{hplot}

