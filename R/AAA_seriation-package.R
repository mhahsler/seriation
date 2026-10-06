#' @keywords internal
#'
#' @section Key functions:
#' - Seriation: [seriate()], [criterion()], [get_order()], [permute()]
#' - Visualization: [pimage()], [bertinplot()], [hmap()], [dissplot()], [VAT()]
#'
#' @section Available seriation methods and criteria:
#' * [A list with the implemented seriation methods](https://michael.hahsler.net/seriation/articles/seriation_methods.html)
#' * [A visual comparison between seriation methods](https://michael.hahsler.net/seriation/articles/comparison.html)
#' * [A list with the implemented seriation criteria](https://michael.hahsler.net/seriation/articles/seriation_criteria.html)
#'
#' @section Quickstart guides:
#' * [How to reorder heatmaps](https://michael.hahsler.net/seriation/articles/heatmaps.html)
#' * [How to reorder correlation matrices](https://michael.hahsler.net/seriation/articles/correlation_matrix.html)
#' * [How to evaluate clusters using dissimilarity plots](https://michael.hahsler.net/seriation/articles/clustering.html)
#'
#' @references Michael Hahsler, Kurt Hornik, and Christian Buchta. 
#' Getting things in order: An introduction to the R package seriation. 
#' Journal of Statistical Software, 25(3):1--34, March 2008. 
#' \doi{10.18637/jss.v025.i03}
#'
#' @importFrom graphics plot text title
#' @importFrom ca ca
#' @importFrom utils head tail
#' @importFrom stats reorder as.dist hclust runif rnorm dist order.dendrogram prcomp
#' @useDynLib seriation, .registration=TRUE
"_PACKAGE"


