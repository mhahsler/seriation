#' @title `r packageDescription("seriation")$Package`: `r packageDescription("seriation")$Title`
#'
#' @description `r packageDescription("seriation")$Description`
#'
#' @section Key functions:
#' - Seriation: [seriate()], [criterion()], [get_order()], [permute()]
#' - Visualization: [pimage()], [bertinplot()], [hmap()], [dissplot()], [VAT()]
#'
#' @author Michael Hahsler
#' @docType package
#' @name seriation
#'
#' @importFrom graphics plot text title
#' @importFrom stats reorder as.dist hclust runif rnorm dist order.dendrogram prcomp
#' @useDynLib seriation, .registration=TRUE
NULL
