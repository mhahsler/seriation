\name{dissplot}
\alias{dissplot}
\alias{plot.reordered_cluster_dissimilarity_matrix}
\alias{print.reordered_cluster_dissimilarity_matrix}
\title{Dissimilarity Plot}
\description{
Visualizes a dissimilarity matrix using seriation and matrix shading
using the method developed by Hahsler and Hornik (2011).
Entries with lower dissimilarities (higher similarity) are plotted darker.
Such a plot can be used to uncover hidden structure in the data.

The plot can also be used to visualize cluster quality (see Ling 1973).
Objects belonging to the same cluster are displayed in consecutive order. The
placement of clusters and the within cluster order is obtained by a seriation
algorithm which tries to place large similarities/small dissimilarities close to
the diagonal. Compact clusters are visible as dark squares (low dissimilarity)
on the diagonal of the plot. Additionally, a Silhouette plot (Rousseeuw 1987)
is added. This visualization is similar to CLUSION (see Strehl and Ghosh 2002),
however, allows for using arbitrary seriating algorithms.
}
\usage{
dissplot(x, labels = NULL, method = "ARSA", 
  control = NULL, options = NULL, \ldots)   
}
\arguments{
  \item{x}{ an object of class \code{dist}.}
  \item{labels}{ \code{NULL} or an integer vector of the same length as
      rows/columns in \code{x} indicating the cluster membership for each
      object in \code{x} as consecutive integers starting with one. The labels
      are used to reorder the matrix.}
  \item{method}{ a list with up to three elements or a single character string.
	  Use a single character string to apply the same algorithm to reorder
	  the clusters (inter cluster seriation) as well as the objects within each
	  cluster (intra cluster seriation).

      If separate algorithms for inter and intra cluster seriation are
      required, \code{method} can be a \code{list} of two named elements
      (\code{inter_cluster} and \code{intra_cluster} each containing the name
      of the respective seriation method. See \code{seriate.dist} for available
      algorithms. 

      Set method to \code{NA} to plot the matrix as is (no or only coarse 
	  seriation). For intra cluster reordering the special
      method \code{"silhouette width"} is available. Objects in clusters are
      then ordered by silhouette width (from silhouette plots).
      If no \code{method} is given, the default method of \code{seriate.dist}
      is used.
      
	  The third list element (named \code{aggregation}) 
	  controls how inter cluster dissimilarities are computed from from
	  the given dissimilarity matrix. The choices are
	   \code{"avg"} (average pairwise dissimilarities; average-link), 
	   \code{"min"} (minimal pairwise dissimilarities; single-link), 
	   \code{"max"} (maximal pairwise dissimilarities; complete-link), and
	   \code{"Hausdorff"} (pairs up each point from one cluster with the most
	   similar point from the other cluster and then uses the largest
	   dissimilarity of paired up points). 
	  }
  \item{control}{ a list of control options passed on to the seriation
      algorithm.  
	  In case of two different seriation algorithms, \code{control} can
      contain a list of two named elements (\code{inter_cluster} 
      and \code{intra_cluster}) containing each a list with the control options
      for the respective algorithm.}
  \item{options}{a list with options for plotting the matrix. The
      list can contain the following elements:

  \describe{
    \item{\code{plot}}{a logical indicating  if a plot should be produced.  if
      \code{FALSE}, the returned object can be plotted later using the function
      \code{plot} which takes as the second argument a list of plotting options
      (see \code{options} below).}
    \item{\code{cluster_labels}}{ a logical indicating whether to display cluster
      labels in the plot.}
    \item{\code{averages}}{ a logical vector of length two. 
	The first element controls the upper triangle and the second element 
	the lower triangle of the plot. \code{FALSE} displays the original
	dissimilarity between objects, \code{TRUE} displays cluster-wise 
	average dissimilarities, and \code{NA} leaves the triangle white
	(default: \code{c(FALSE, TRUE)}, i.e., the lower triangle displays
	averages)}
    \item{\code{lines}}{ a logical indicating  whether to draw lines to separate
      clusters.} 
    \item{\code{flip}}{ a logical indicating if the clusters are displayed 
  on the diagonal from
  north-west to south-east (\code{FALSE}; default) or from 
  north-east to south-west (\code{TRUE}).} 
    \item{\code{silhouettes}}{ a logical indicating  whether to include a
      silhouette plot (see Rousseeuw, 1987).} 
    \item{\code{threshold}}{ a numeric. If used, only plot distances below the
      threshold are displayed. Consider also using \code{zlim} for this purpose.} 
    \item{\code{col}}{ colors used for the image plot.}
    \item{\code{key}}{ a logical indicating  whether to place a color key
      below the plot. }
    \item{\code{zlim}}{ range of values to display (defaults to range \code{x}).}
    \item{\code{axes}}{ \code{"auto"} (default; enabled for less than 
        25 objects), \code{"y"} or \code{"none"}.}
    \item{\code{main}}{ title for the plot.}
    \item{\code{newpage}}{ a logical indicating  whether to start plot on a new
      page (see \code{grid.newpage} in package \pkg{grid}). }
    \item{\code{pop}}{a logical indicating  whether to pop the created viewports
      (see package \pkg{grid})?  }
    \item{\code{gp}, \code{gp_lines}, \code{gp_labels}}{ 
      objects of class \code{gpar}       
      containing graphical
      parameters (see \code{gpar} in package \pkg{grid}). }
      }
  }
  \item{\ldots}{further arguments are added to \code{options}.}
}
\value{
An invisible object of class \code{cluster_proximity_matrix} with the following
elements:
\item{order}{ \code{NULL} or integer vector giving the order used to plot
    \code{x}.}
\item{cluster_order}{ \code{NULL} or integer vector giving the order 
    of the clusters as plotted.}
\item{method}{ vector of character strings indicating the seriation methods 
    used for plotting  \code{x}.}
\item{k}{ \code{NULL} or integer scalar giving the number of clusters
    generated.}
\item{description}{ a \code{data.frame}  containing information (label, size,
    average intra-cluster dissimilarity and the average silhouette) for the
    clusters as displayed in the plot (from top/left to bottom/right).}

This object can be used for plotting via
\code{plot(x, options = NULL, ...)}, where \code{x} is the
object and \code{options} contains a list with plotting options (see above).
}
\references{
Hahsler, M. and Hornik, K. (2011):
 Dissimilarity plots: A visual exploration tool for partitional clustering. \emph{Journal of Computational and Graphical Statistics,} \bold{10}(2):335--354.

Ling, R.F. (1973):
A computer generated aid for cluster analysis. 
\emph{Communications of the ACM,} 
\bold{16}(6), 355--361.

Rousseeuw, P.J. (1987): 
Silhouettes: A graphical aid to the interpretation and 
validation of cluster analysis. 
\emph{Journal of Computational and Applied Mathematics,}
\bold{20}(1), 53--65.

Strehl, A. and Ghosh, J. (2003):
Relationship-based clustering and 
visualization for high-dimensional data mining. 
\emph{INFORMS Journal on Computing,} 
\bold{15}(2), 208--230.

}
\seealso{
  \code{\link[stats]{dist}},
  \code{\link{seriate}},
  \code{\link{pimage}} and
  \code{\link{hmap}}.
}
\author{Michael Hahsler}
\examples{
data("iris")
d <- dist(iris[-5])

## plot original matrix
res <- dissplot(d, method = NA)

## plot reordered matrix using the nearest insertion algorithm (from tsp)
res <- dissplot(d, method = "TSP",
    options = list(main = "Seriation (TSP)"))

## cluster with pam (we know iris has 3 clusters)
library("cluster")
l <- pam(d, 3, cluster.only = TRUE)

## we use a grid layout to place several plots on a page
library("grid")
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow = 2, ncol = 2), 
    gp = gpar(fontsize = 8)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

## visualize the clustering (using ARSA between clusters and MDS within)
res <- dissplot(d, l, method = list(inter = "ARSA", intra = "MDS"),  
    options = list(main = "PAM + Seriation - standard", 
    newpage = FALSE))

popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))

## more visualization options. Note that we reuse the reordered object res!
## color: use 10 shades red-blue
plot(res, options = list(main = "PAM + Seriation", 
    col= bluered(10, bias=.5), newpage = FALSE))

popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))

## threshold (using zlim) and cubic scale to highlight differences
plot(res, options = list(main = "PAM + Seriation - threshold", 
    zlim = c(0, 1.5), col = greys(100, power = 2), newpage = FALSE))

popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))

## use custom (logistic) scale
plot(res, options = list(main = "PAM + Seriation - logistic scale", 
    col= hcl(c = 0, l = (plogis(seq(10, 0, length=100), 
	location = 2, scale = 1/2, log = FALSE))*100), 
	newpage = FALSE))

popViewport(2)

## the reordered_cluster_dissimilarity_matrix object
res 
names(res)
}
\keyword{hplot}
\keyword{cluster}
