\name{reorder.hclust}
\alias{reorder}
\alias{reorder.hclust}
\title{Reorder Dendrograms using Optimal Leaf Ordering}
\description{
Reorder method for dendrograms for optimal leaf ordering.
}
\usage{
\method{reorder}{hclust}(x, dist, method = "OLO", ...)
}
\arguments{
  \item{x}{an object of class \code{hclust}.}
  \item{dist}{an object of class \code{dist} with dissimilarities between the
    objects in \code{x}.}
  \item{method}{ a character string with the name of the used measure
    Available are: 
      \code{"OLO"} (optimal leaf ordering; Bar-Joseph et al., 2001) and
      \code{"GW"} (Gruvaeus and Wainer, 1972).}
  \item{...}{ further arguments are currently ignored.}
}
\details{
Minimizes the distance between neighboring objects
(leaf nodes) in the dendrogram by flipping the order of subtrees. 
The algorithm by Gruvaeus and Wainer is implemented in
package \pkg{gclus} (Hurley 2004).
}
\value{
A reordered \code{hclust} object.
}
\references{ 
Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola. (2001): Fast
Optimal Leaf Ordering for Hierarchical Clustering.
\emph{Bioinformatics,} \bold{17}(1), 22--29.
  
Gruvaeus, G. and Wainer, H. (1972): Two Additions to Hierarchical Cluster
Analysis,
\emph{British Journal of Mathematical and Statistical Psychology,}
\bold{25}, 200--206.
  
Hurley, Catherine B. (2004): Clustering Visualizations of Multidimensional
Data.
\emph{Journal of Computational and Graphical Statistics,}
\bold{13}(4), 788--806.
}
\seealso{
\code{\link[gclus]{reorder.hclust}}
}
\author{Michael Hahsler}
\examples{
## cluster European cities by distance
data("eurodist")
d <- as.dist(eurodist)
hc <- hclust(eurodist)
  
## plot original dendrogram and the reordered dendrograms  
plot(hc)  
plot(reorder(hc, d, method = "GW"))  
plot(reorder(hc, d, method = "OLO"))  
}
\keyword{optimize}
\keyword{cluster}
