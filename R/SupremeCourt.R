#' Voting Patterns in the Second Rehnquist U.S. Supreme Court
#'
#' Contains a (a subset of the) decisions for the stable 8-yr
#' period 1995-2002 of the second Rehnquist Supreme Court.
#' Decisions are aggregated to
#' the joint probability for disagreement between judges.
#'
#' @name SupremeCourt
#' @aliases SupremeCourt
#' @docType data
#' @family data
#' @format
#'   A square, symmetric 9-by-9 matrix with the joint probability for disagreement.
#' @references
#'     Sirovich, L. (2003). A pattern analysis of the second Rehnquist
#'     U.S. Supreme Court. _Proceedings of the National Academy of Sciences of the United
#'     States of America,_ **100**, 7432-7437. \doi{10.1073/pnas.1132164100}
#' @author Michael Hahsler
#' @examples
#' data("SupremeCourt")
#'
#' # a matrix with joint probability of disagreement
#' SupremeCourt
#'
#' # show judges in original alphabetical order
#' d <- as.dist(SupremeCourt)
#' pimage(d, diag = TRUE, upper = TRUE)
#'
#' # reorder judges using seriation based on similar decisions
#' o <- seriate(d)
#' o
#'
#' pimage(d, o, diag = TRUE, upper = TRUE)
#'
#' # Use optimal leaf ordering (hierarchical clustering with reordering)
#' # which uses a dendrogram
#' o <- seriate(d, method = "OLO")
#' o
#'
#' plot(o[[1]])
#'
#' # Use multi-dimensional scaling and show the configuration
#' o <- seriate(d, method = "sammon")
#' o
#'
#' pimage(d, o, diag = TRUE, upper = TRUE)
#' plot_config(o[[1]])
#' @keywords datasets
NULL
