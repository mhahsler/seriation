\name{hmap}
\alias{hmap}
\title{Plot Heat Map Reordered Using Seriation}
\description{
Provides heat maps reordered using several different seriation methods.
This includes dendrogram based reordering with optimal leaf order and matrix seriation based heat maps.
}
\usage{
hmap(x, distfun = dist, method = "OLO", control = NULL, zlim = NULL, ...)
}
\arguments{
  \item{x}{a matrix or a dissimilarity matrix of class dist. If a dissimilarity
    matrix is used, then the \code{distfun} is ignored. }
  \item{distfun}{function used to compute the distance (dissimilarity) between
      both rows and columns (default: \code{dist}). }
  \item{method}{a character strings indicating the used seriation algorithm 
      (see \code{seriate.dist}). If the method results in a dendrogram then
      \code{heatmap} in \pkg{stats} is used to show the dendrograms, otherwise
      reordered distance matrices are shown instead.}
  \item{control}{a list of control options passed on to the seriation 
      algorithm specified in \code{method}.} 
  \item{zlim}{range of values to display (defaults to the range 
      of \code{x}).}
  \item{\dots}{further arguments.} 
}  
\details{
  
  For dendrogram based heat maps the 
      arguments are passed on to \code{heatmap.2} in
      \pkg{gplots}. See for example \code{margins} and \code{col}.
      The following arguments for \code{heatmap.2} cannot be used:
      \code{Rowv, Colv, hclustfun, reorderfun}.
      
  For seriation-based heat maps further arguments include:
       \describe{
          \item{\code{gp}}{an object of class \code{gpar} containing graphical
                parameters (see \code{gpar} in package \pkg{grid}).}
          \item{\code{newpage}}{a logical indicating  whether to start 
                plot on a new
                page (see \code{gpar} in package \pkg{grid}).} 
          \item{\code{prop}}{a logical indicating whether the height and width of
                \code{x} should be plotted proportional to its dimensions.}
          \item{\code{showdist}}{Display seriated dissimilarity matrices?                      Values are \code{"none"}, \code{"both"}, 
                \code{"rows"} or \code{"columns"}. } 
         \item{\code{key}}{logical; show a colorkey?.} 
         \item{key.lab}{ string plotted next to the color key.}
         \item{\code{axes}}{ one of \code{"auto"} (default; show axis 
                labels if there are
                less than 25 labels), \code{"x"}, \code{"y"}, 
                \code{"both"} and \code{"none"}.} 
         \item{\code{margins}}{ bottom and right-hand-side margins are 
                calculated automatically or can be specifies as a vector 
                of two numbers (in lines).} 
         \item{\code{zlim}}{ range of values displayed.} 
         \item{\code{col}, \code{col_dist}}{ color palettes used.} 
      }      

For \code{dendrogram = TRUE}, \code{seriate.hclust} 
with the default method \code{"optimal"} is used for arranging the
dendrograms and \code{x}. \code{heatmap} is used for plotting.

For \code{dendrogram = FALSE}, \code{seriate.dist} 
with the default method \code{"tsp"}
(a traveling salesperson solver) for arranging \code{x} is used.
\pkg{grid} code implemented in this package 
is used to produce the plot.

Note that unlike the default behavior of \code{heatmap}, 
scaling is not automatically applied. The data
have to be scaled before using \code{hmap}.
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
\code{\link[gplots]{heatmap.2}} in \pkg{gplots}.}
\author{Michael Hahsler}
\examples{
data("Wood")

## default heatmap does Euclidean distance, hierarchical clustering with
## average-link and optimal leaf ordering
hmap(Wood)  
  
## heatmap with correlation-based distance, green-red color (greenred is 
## predefined) and optimal leaf ordering and no row label
dist_cor <- function(x) as.dist(1-cor(t(x)))
  
hmap(Wood, method="OLO", distfun = dist_cor, col=greenred(100), labRow=FALSE)
  
## order-based heatmap
hmap(Wood, method="MDS_angle", distfun = dist_cor, col=greenred(100))  

## order-based with dissimilarity matrices
hmap(Wood, method="MDS_angle", distfun = dist_cor, showdist = "both", 
  col=greenred(100))  
}
\keyword{hplot}
