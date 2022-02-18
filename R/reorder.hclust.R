#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



#' Reorder Dendrograms using Optimal Leaf Ordering
#'
#' Reorder method for dendrograms for optimal leaf ordering.
#'
#' Minimizes the distance between neighboring objects (leaf nodes) in the
#' dendrogram by flipping the order of subtrees. The algorithm by Gruvaeus and
#' Wainer is implemented in package \pkg{gclus} (Hurley 2004).
#'
#' @aliases reorder reorder.hclust
#' @param x an object of class \code{hclust}.
#' @param dist an object of class \code{dist} with dissimilarities between the
#' objects in \code{x}.
#' @param method a character string with the name of the used measure.
#' Available are:
#'   - \code{"OLO"} (optimal leaf ordering; Bar-Joseph et al., 2001) implemented in this package and
#'   - \code{"GW"} (Gruvaeus and Wainer, 1972) from package \pkg{gclus}.
#' @param ...  further arguments are currently ignored.
#' @return A reordered \code{hclust} object.
#' @author Michael Hahsler
#' @seealso [gclus::reorder.hclust()]
#' @references Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola.
#' (2001): Fast Optimal Leaf Ordering for Hierarchical Clustering.
#' \emph{Bioinformatics,} \bold{17}(1), 22--29.
#'
#' Gruvaeus, G. and Wainer, H. (1972): Two Additions to Hierarchical Cluster
#' Analysis, \emph{British Journal of Mathematical and Statistical Psychology,}
#' \bold{25}, 200--206.
#'
#' Hurley, Catherine B. (2004): Clustering Visualizations of Multidimensional
#' Data. \emph{Journal of Computational and Graphical Statistics,}
#' \bold{13}(4), 788--806.
#' @keywords optimize cluster
#' @examples
#' ## cluster European cities by distance
#' data("eurodist")
#' d <- as.dist(eurodist)
#' hc <- hclust(eurodist)
#'
#' ## plot original dendrogram and the reordered dendrograms
#' plot(hc)
#' plot(reorder(hc, d, method = "GW"))
#' plot(reorder(hc, d, method = "OLO"))
#' @export
reorder.hclust <- function(x, dist, method = "OLO", ...) {
  method <- match.arg(tolower(method), choices = c("olo", "gw"))

  ## no reordering for less than 3 objects!
  if (length(x$order) < 3)
    return(x)

  switch(method,
    olo = .seriate_optimal(x, dist),
    gw  = .seriate_gruvaeus(x, dist))
}

## wrapper for reorder.hclust in gclus
.seriate_gruvaeus <- function(hclust, dist)
  gclus::reorder.hclust(hclust, dist)

## wrapper to the optimal leaf ordering algorithm
##
## ceeboo 2005
.seriate_optimal <- function(hclust, dist) {
  ## check hclust
  merge <- hclust$merge
  if (!is.matrix(merge))
    stop("Component 'merge' of argument 'hclust' must be a matrix.")
  if (length(dim(merge)) != 2)
    stop("Component 'merge' of argument 'hclust' is invalid.")
  if (dim(merge)[1] != attr(dist, "Size") - 1)
    stop("Argument 'dist' and component 'merge' of argument 'hclust' do not conform.")
  mode(merge) <- "integer"

  obj <- .Call("order_optimal", dist, merge)

  names(obj) <- c("merge", "order", "length")
  ##names(obj$order) <- attr(dist,"Labels")
  hclust$merge <- obj$merge
  hclust$order <- obj$order

  hclust
}
