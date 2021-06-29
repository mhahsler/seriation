\name{hmap}
\alias{hmap}
\alias{gghmap}
\alias{heatmap}
\title{Plot Heat Map Reordered Using Seriation}
\description{
Provides heatmaps reordered using several different seriation methods.
This includes dendrogram based reordering with optimal leaf order and matrix seriation based heat maps.
}
\usage{
hmap(x, distfun = dist, method = "OLO", control = NULL,
  scale = c("none", "row", "column"), showDend = TRUE, col = NULL,
  row_labels = NULL, col_labels = NULL, ...)

gghmap(x, distfun = dist, method = "OLO", control = NULL,
  scale = c("none", "row", "column"), prop = FALSE, ...)
}
\arguments{
  \item{x}{a matrix or a dissimilarity matrix of class dist. If a dissimilarity
    matrix is used, then the \code{distfun} is ignored. }
  \item{distfun}{function used to compute the distance (dissimilarity) between
      both rows and columns (default: \code{dist}). For \code{gghmap}, this parameter is
      passed on in \code{control}. }
  \item{method}{a character strings indicating the used seriation algorithm
      (see \code{seriate.dist}). If the method results in a dendrogram then
      \code{heatmap} in \pkg{stats} is used to show the dendrograms, otherwise
      reordered distance matrices are shown instead.}
  \item{control}{a list of control options passed on to the seriation
      algorithm specified in \code{method}.}
  \item{scale}{character indicating if the values should be centered and scaled in either
    the row direction or the column direction, or none. Default is none.}
  \item{showDend}{ Show dendrograms in the margin? }
  \item{col}{a list of colors used.}
  \item{row_labels, col_labels}{ a logical indicating if row and column labels in \code{x} should be displayed.
    If \code{NULL} then labels are displayed if the \code{x} contains the appropriate dimname and the number of labels is 25 or less. A character vector of the appropriate length with labels can also be supplied. }
  \item{prop}{logical; change the aspect ratio so cells in the image have a equal width and height.}
  \item{\dots}{further arguments passed on to \code{\link{heatmap}} in package \pkg{stats}.}
}
\details{

  For dendrogram based heat maps, the
      arguments are passed on to \code{heatmap} in
      \pkg{stats}.
      The following arguments for\code{heatmap} cannot be used:
      \code{margins, Rowv, Colv, hclustfun, reorderfun}.

  For seriation-based heat maps further arguments include:
       \describe{
          \item{\code{gp}}{an object of class \code{gpar} containing graphical
                parameters (see \code{gpar} in package \pkg{grid}).}
          \item{\code{newpage}}{a logical indicating  whether to start
                plot on a new
                page (see \code{gpar} in package \pkg{grid}).}
          \item{\code{prop}}{a logical indicating whether the height and width of
                \code{x} should be plotted proportional to its dimensions.}
          \item{\code{showdist}}{Display seriated dissimilarity matrices?
              Values are \code{"none"}, \code{"both"},
                \code{"rows"} or \code{"columns"}. }
         \item{\code{key}}{logical; show a colorkey?.}
         \item{key.lab}{ string plotted next to the color key.}
#         \item{\code{}}{ one of \code{"auto"} (default; show axis
#                labels if there are
#                less than 25 labels), \code{"x"}, \code{"y"},
#                \code{"both"} and \code{"none"}.}
         \item{\code{margins}}{ bottom and right-hand-side margins are
                calculated automatically or can be specifies as a vector
                of two numbers (in lines).}
         \item{\code{zlim}}{ range of values displayed.}
         \item{\code{col}, \code{col_dist}}{ color palettes used.}
      }
}
\value{
An invisible list with elements:
\item{rowInd, colInd}{index permutation vectors.}
\item{reorder_method}{name of the method used to reorder the matrix.}

The list may contain additional elements (dendrograms, colors, etc).

    }
\seealso{
\code{\link{seriate}},
\code{\link{pimage}},
\code{\link{dissplot}},
\code{\link[stats]{heatmap}} in \pkg{stats}.}
\author{Michael Hahsler}
\examples{
data("Wood")

# regular heatmap from package stats
heatmap(Wood, main = "Wood (standard heatmap)")

# Default heatmap does Euclidean distance, hierarchical clustering with
# average-link and optimal leaf ordering
hmap(Wood, main = "Wood (opt. leaf ordering)")

# Heatmap shown without dendrograms (used pimage)
hmap(Wood, main = "Wood (opt. leaf ordering)",  showDend = FALSE)

# Heatmap with correlation-based distance, green-red color (greenred is
# predefined) and optimal leaf ordering and no row label
dist_cor <- function(x) as.dist(1 - cor(t(x)))
hmap(Wood, distfun = dist_cor, col = greenred(100), row_labels = FALSE,
  main = "Wood (reorded by corr. between obs.)")

# Heatmap with order based on the angle in two-dimensional MDS space.
hmap(Wood, method = "MDS_angle", col = greenred(100), row_labels = FALSE,
  main = "Wood (reorderd using ange in MDS space)")

# Heatmap for distances
d <- dist(Wood)
hmap(d, method = "OLO", main = "Wood (Euclidean distances)")

# order-based with dissimilarity matrices
hmap(Wood, method = "MDS_angle", showdist = "both",
  col = greenred(100), col_dist = greens(100),
  keylab = "norm. Expression", main = "Wood (reporderd with distances)")

# Manually seriate and plot as pimage.
o <- seriate(Wood, method = "heatmap", control = list(dist_fun = dist, seriation_method = "OLO"))
o

pimage(Wood, o, prop = FALSE)

# Note: method heatmap calculates reorderd hclust objects which can be used for many heatmap
# implementations.
heatmap(Wood, Rowv = as.dendrogram(o[[1]]), Colv = as.dendrogram(o[[2]]))

# ggplot 2 version does not support dendrograms in the margin (for now)
if (require("ggplot2")) {
  library("ggplot2")

  gghmap(Wood) + labs(title = "Wood", subtitle = "Optimal leaf ordering")

  # More parameters (see ? ggpimage): reverse column order and flip axes, make a proportional plot
  gghmap(Wood, reverse_columns = TRUE) +
    labs(title = "Wood", subtitle = "Optimal leaf ordering")

  gghmap(Wood, flip_axes = TRUE) +
    labs(title = "Wood", subtitle = "Optimal leaf ordering")

  gghmap(Wood, flip_axes = TRUE, prop = TRUE) +
    labs(title = "Wood", subtitle = "Optimal leaf ordering")

  gghmap(Wood, distfun = dist_cor) +
    labs(title = "Wood", subtitle = "Reorded by correlation between observations") +
    scale_fill_gradient2(low = "darkgreen", high = "red")

  gghmap(d, prop = TRUE) +
    labs(title = "Wood", subtitle = "Euclidean distances, reordered)")

  # Note: the ggplot2-based version cannot show distance matrices in the same plot.

  # Manually seriate and plot as pimage.
  o <- seriate(Wood, method = "heatmap", control = list(dist_fun = dist,
    seriation_method = "OLO"))
  o

  ggpimage(Wood, o, prop = FALSE)
}

}
\keyword{hplot}
