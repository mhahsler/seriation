#' Extracting Order Information from a Permutation Object
#'
#' Method to get the order information from an object of class
#' [ser_permutation] or [ser_permutation_vector]. Order information
#' can be extracted as a permutation vector, a vector containing each
#' object's rank or a permutation matrix.
#'
#' `get_order()` returns the permutation as an integer vector which arranges the
#' objects in the seriation order. That is, a vector with the index of the first,
#' second, \eqn{..., n}-th object in the order defined by the permutation.
#' These permutation vectors can directly be
#' used to reorder objects using subsetting with `"["`.  \emph{Note:} In
#' \pkg{seriation} we usually use these order-based permutation vectors.
#' **Note on names:** While R's [order()] returns an unnamed vector,
#' `get_order()` returns names (if available). The names are the object label
#' corresponding to the index at that position.
#' Therefore, the names in the order are in the order after
#' the permutation.
#'
#' `get_rank()` returns the seriation as an integer vector containing the
#' rank/position for each objects after the permutation is applied.
#' That is, a vector with the position of the first, second,
#' \eqn{..., n}-th object after permutation.  Note: Use
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
#' ## create a random ser_permutation_vector
#' ## Note that ser_permutation_vector is a single permutation vector
#' x <- structure(1:10, names = paste0("X", 1:10))
#' o <- sample(x)
#' o
#'
#' p <- ser_permutation_vector(o)
#' p
#'
#' get_order(p)
#' get_rank(p)
#' get_permutation_matrix(p)
#'
#' ## reorder objects using subsetting, the provided permute function or by
#' ## multiplying the with the permutation matrix. We use here
#' x[get_order(p)]
#' permute(x, p)
#' drop(get_permutation_matrix(p) %*%  x)
#'
#' ## ser_permutation contains one permutation vector for each dimension
#' p2 <- ser_permutation(p, sample(5))
#' p2
#'
#' get_order(p2, dim = 2)
#' get_rank(p2, dim = 2)
#' get_permutation_matrix(p2, dim = 2)
#' @export
get_order <- function(x, ...)
  UseMethod("get_order")

#' @export
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

#' @rdname get_order
#' @export
get_order.numeric <- function(x, ...) {
  structure(order(x), names = names(x))
}

## returns for each object its rank (rank of first, second, etc. object)
#' @rdname get_order
#' @export
get_rank <- function(x, ...) {
  o <- get_order(x, ...)
  r <- order(o)
  names(r) <- names(o)[r]

  r
}


#' @rdname get_order
#' @export
get_permutation_matrix <- function(x, ...)
  permutation_vector2matrix(get_order(x, ...))
