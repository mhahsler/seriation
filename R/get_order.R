#' Extracting Order Information from a Permutation Object
#'
#' Method to get the order information from an object of class
#' [ser_permutation] or [ser_permutation_vector]. Order information
#' can be extracted as an integer permutation vector, a vector containing the
#' object ranks or a permutation matrix.
#'
#' `get_order()` returns the seriation as an integer vector containing the
#' order of the objects after permutation. That is, the index of the first,
#' second, \eqn{..., n}-th object. These permutation vectors can directly be
#' used to reorder objects using subsetting with `"["`.  \emph{Note:} In
#' \pkg{seriation} we usually use these order-based permutation vectors.
#'
#' `get_rank()` returns the seriation as an integer vector containing the
#' rank/position for each objects in the permutation. That is, position of the
#' first, second, \eqn{..., n}-th object after permutation.  Note: Use
#' `order()` to convert ranks back to an order.
#'
#' `get_permutation_matrix()` returns a \eqn{n \times n}{n x n} permutation
#' matrix.
#'
#' @family permutation
#'
#' @param x an object of class [ser_permutation] or
#' [ser_permutation_vector].
#' @param dim order information for which dimension should be returned?
#' @param ... further arguments are ignored for `get_order()`.  For
#' `get_rank()` and for `get_permutation_matrix()` the additional
#' arguments are passed on to `get_order()` (e.g., as `dim`).
#' @return Returns an integer permutation vector/a permutation matrix.
#'
#' @author Michael Hahsler
#' @keywords manip
#' @examples
#' ## ser_permutation_vector is a single permutation vector
#' o <- ser_permutation_vector(sample(10))
#' o
#'
#' get_order(o)
#' get_rank(o)
#' get_permutation_matrix(o)
#'
#' ## ser_permutation contains one permutation vector for each dimension
#' o2 <- ser_permutation(o, sample(5))
#' o2
#'
#' get_order(o2, dim = 2)
#' get_rank(o2, dim = 2)
#' get_permutation_matrix(o2, dim = 2)
#' @export
get_order <- function(x, ...)
  UseMethod("get_order")

get_order.default <- function(x, ...)
  stop(gettextf("No permutation accessor implemented for class '%s'. ",
    class(x)))

#' @rdname get_order
#' @export
get_order.ser_permutation_vector <- function(x, ...)
  NextMethod()


#' @rdname get_order
#' @export
get_order.ser_permutation <-
  function(x, dim = 1, ...)
    get_order(x[[dim]])


#' @rdname get_order
#' @export
get_order.hclust <- function(x, ...)
  structure(.Data = x$order, names = x$labels[x$order])


#' @rdname get_order
#' @export
get_order.dendrogram <- function(x, ...)
  order.dendrogram(x)

#' @rdname get_order
#' @export
get_order.integer <- function(x, ...) {
  if (.is_identity_permutation(x))
    stop("Cannot get order vector from symbolic identity permutation (undefined length).")
  structure(as.integer(x), names = names(x))
}


## returns for each object its rank (rank of first, second, etc. object)
#' @rdname get_order
#' @export
get_rank <- function(x, ...)
  order(get_order(x, ...))

#' @rdname get_order
#' @export
get_permutation_matrix <- function(x, ...)
  permutation_vector2matrix(get_order(x, ...))
