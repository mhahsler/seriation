#' @title `r packageDescription("seriation")$Package`: `r packageDescription("seriation")$Title`
#'
#' @description Infrastructure for ordering objects with an implementation of several seriation/sequencing/ordination techniques to reorder matrices, dissimilarity matrices, and dendrograms. Also provides (optimally) reordered heatmaps, color images and clustering visualizations like dissimilarity plots, and visual assessment of cluster tendency plots (VAT and iVAT).
#'
#' @references Michael Hahsler, Kurt Hornik, and Christian Buchta. Getting things in order: An introduction to the R package seriation. Journal of Statistical Software, 25(3):1--34, March 2008. \doi{10.18637/jss.v025.i03}
#'
#'
#' @section Key functions:
#' - Seriation: [seriate()], [criterion()], [get_order()], [permute()]
#' - Visualization: [pimage()], [bertinplot()], [hmap()], [dissplot()], [VAT()]
#'
#' @author Michael Hahsler
#' @docType package
#' @name seriation-package
#'
#' @importFrom graphics plot text title
#' @importFrom stats reorder as.dist hclust runif rnorm dist order.dendrogram prcomp
#' @useDynLib seriation, .registration=TRUE
NULL
