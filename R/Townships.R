#' Bertin's Characteristics of Townships
#'
#' This data contains nine characteristics for 16 townships. The data
#' set was used by Bertin (1981) to illustrate that the conciseness
#' of presentation can be improved by seriating the rows and columns.
#'
#' @name Townships
#' @aliases Townships
#' @family data
#' @docType data
#' @format
#'   A matrix with 16 0-1 variables (columns) indicating the presence
#'   (`1`) or absence (`0`) of characteristics of townships
#'   (rows).
#' @references
#' Bertin, J. (1981): _Graphics and Graphic Information Processing_. Berlin, Walter de Gruyter.
#' @author Michael Hahsler
#' @examples
#' data("Townships")
#'
#' ## original data
#' pimage(Townships)
#' criterion(Townships)
#'
#' ## seriated data
#' order <- seriate(Townships, method = "BEA", control = list(rep = 5))
#' pimage(Townships, order)
#' criterion(Townships, order)
#' @keywords datasets
NULL
