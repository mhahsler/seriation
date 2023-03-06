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


# helper
ndim <- function(x)
  length(dim(x))

#' Permute the Order in Various Objects
#'
#' Provides the generic function and methods for permuting the order of various
#' objects including vectors, lists, dendrograms (also \code{hclust} objects),
#' the order of observations in a \code{dist} object, the rows and columns of a
#' matrix or data.frame, and all dimensions of an array given a suitable
#' [ser_permutation] object.
#'
#' The permutation vectors in [ser_permutation] are suitable if the number
#' of permutation vectors matches the number of dimensions of \code{x} and if
#' the length of each permutation vector has the same length as the
#' corresponding dimension of \code{x}.
#'
#' For 1-dimensional/1-mode data (list, vector, \code{dist}), \code{order} can
#' also be a single permutation vector of class [ser_permutation_vector]
#' or data which can be automatically coerced to this class (e.g. a numeric
#' vector).
#'
#' For \code{dendrogram} and \code{hclust}, subtrees are rotated to represent
#' the order best possible. If the order is not achieved perfectly then the
#' user is warned. This behavior can be changed with the extra parameter
#' \code{incompatible} which can take the values \code{"warn"} (default),
#' \code{"stop"} or \code{"ignore"}.
#'
#' @family permutation
#'
#' @param x an object (a list, a vector, a \code{dist} object, a matrix, an
#' array or any other object which provides \code{dim} and standard subsetting
#' with \code{"["}).
#' @param order an object of class [ser_permutation] which contains
#' suitable permutation vectors for \code{x}. Alternatively, a character string with the
#' name of a seriation method appropriate for `x` can be specified (see [seriate()]).
#' This will perform seriation and permute `x`.
#' @param margin specifies the dimensions to be permuted as a vector with dimension indices.
#' If `NULL`, \code{order} needs to contain a permutation for all dimensions.
#' If a single margin is specified, then \code{order} can also contain
#' a single permutation vector.
#' \code{margin} are ignored.
#' @param ...  if `order` is the name of a seriation method, then additional arguments are
#' passed on to [seriate()].
#' @returns A permuted object of the same class as `x`.
#' @author Michael Hahsler
#' @keywords manip
#' @examples
#' # List data types for permute
#' methods("permute")
#'
#' # Permute matrix
#' m <- matrix(rnorm(10), 5, 2, dimnames = list(1:5, LETTERS[1:2]))
#' m
#'
#' # Permute rows and columns
#' o <- ser_permutation(5:1, 2:1)
#' o
#'
#' permute(m, o)
#'
#' ## permute only columns
#' permute(m, o, margin = 2)
#'
#' ## permute using PCA seriation
#' permute(m, "PCA")
#'
#' ## permute only rows using PCA
#' permute(m, o, margin = 1)
#'
#' # Permute data.frames
#' df <- as.data.frame(m)
#' permute(df, o)
#'
#' # Permute objects in a dist object
#' d <- dist(m)
#' d
#'
#' permute(d, c(3, 2, 1, 4, 5))
#'
#' permute(d, "Spectral")
#'
#' # Permute a list
#' l <- list(a = 1:5, b = letters[1:3], c = 0)
#' l
#'
#' permute(l, c(2, 3, 1))
#'
#' # Permute a dendrogram
#' hc <- hclust(d)
#' plot(hc)
#' plot(permute(hc, 5:1))
#' @export
permute <- function(x, order, ...)
  UseMethod("permute")

#' @export
permute.default <- function(x, order, ...)
  .permute_kd(x, order, ...)

#' @rdname permute
#' @export
permute.array <- function(x, order, margin = NULL, ...)
  .permute_kd(x, order, margin = margin, ...)

#' @rdname permute
#' @export
permute.matrix <- function(x, order, margin = NULL, ...)
  .permute_kd(x, order, margin = margin, ...)

#' @rdname permute
#' @export
permute.data.frame <- function(x, order, margin = NULL, ...)
  .permute_kd(x, order, margin = margin, ...)

#' @rdname permute
#' @export
permute.table <- function(x, order, margin = NULL, ...)
  .permute_kd(x, order, margin = margin, ...)

#' @rdname permute
#' @export
permute.numeric <- function(x, order, ...)
  .permute_1d(x, order, ...)

#' @rdname permute
#' @export
permute.character <- function(x, order, ...)
  .permute_1d(x, order, ...)

#' @rdname permute
#' @export
permute.list <- function(x, order, ...)
  .permute_1d(x, order, ...)

# special cases
#' @rdname permute
#' @export
permute.dist <- function(x, order, ...) {
  if (is.character(order))
    order <- seriate(x, method = order, ...)

  .nodots(...)

  if (!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  if (.is_identity_permutation(order[[1]]))
    return(x)

  .check_dist_perm(x, order)

  .rearrange_dist(x, get_order(order, 1))
}

#' @rdname permute
#' @export
permute.dendrogram <- function(x, order, ...) {
  .nodots(...)

  if (length(get_order(order)) != stats::nobs(x))
    stop("Length of order and number of leaves in dendrogram do not agree!")


  # modeled after rotate in dendextend. Copied here to reduce the heavy dependency count of dendextend.
  #  x <- dendextend::rotate(x, order = match(get_order(order), get_order(x)))
  rot <- function (x, order, ...)
  {
    if (missing(order)) {
      warning("'order' parameter is missing, returning the tree as it was.")
      return(x)
    }
    labels_x <- labels(x)
    order_x <- order.dendrogram(x)
    number_of_leaves <- length(order_x)
    if (!is.numeric(order)) {
      order <- as.character(order)
      if (length(intersect(order, labels_x)) != number_of_leaves) {
        stop(
          "'order' is neither numeric nor a vector with ALL of the labels (in the order you want them to be)"
        )
      }
      order <- match(order, labels_x)
    }
    weights <- seq_len(number_of_leaves)
    weights_for_order <- numeric(number_of_leaves)
    weights_for_order[order_x[order]] <- weights
    reorder(x, weights_for_order, mean, ...)
  }

  x <- rot(x, order = match(get_order(order), get_order(x)))

  if (any(get_order(x) != get_order(order)))
    warning("Dendrogram cannot be perfectly reordered! Using best approximation.")

  x
}

#' @rdname permute
#' @export
permute.hclust <- function(x, order, ...) {
  nd <- stats::as.hclust(permute(stats::as.dendrogram(x), order, ...))
  x$merge <- nd$merge
  x$height <- nd$height
  x$order <- nd$order

  x
}

# helper
.check_dist_perm <- function(x, order) {
  if (length(order) != 1L)
    stop("dimensions do not match")

  if (attr(x, "Size") != length(get_order(order, 1)))
    stop("some permutation vectors do not fit dimension of data")

  # check dist
  if (isTRUE(attr(x, "Diag")) || isTRUE(attr(x, "Upper")))
    stop("'dist' with diagonal or upper triangle matrix not implemented")
}

.check_matrix_perm <- function(x, order) {
  if (ndim(x) != length(order))
    stop("dimensions do not match")
  if (any(dim(x) != sapply(order, length)))
    stop("some permutation vectors do not fit dimension of data")
}

.permute_kd <- function(x, order, margin = NULL, ...) {

  # DEPRECATED: Compatibility with old permutation for data.frame
  if (is.data.frame(x) && is.null(margin) && length(order) == 1 && !is.character(order)) {
    message(
      "permute for data.frames with a single seriation order is now deprecated. Specify the margin as follows: 'permute(x, order, margin = 1)'"
    )
    margin <- 1
  }

  if (is.null(margin))
    margin <- seq(ndim(x))
  else {
    margin <- as.integer(margin)
    if (!all(margin %in% seq(ndim(x))))
      stop("all margins need to specify a valid dimension in x")
  }


  if (is.character(order))
    order <- seriate(x, method = order, margin = margin, ...)

  #.nodots(...)

  if (!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  if (length(order) != ndim(x) && length(order) != length(margin))
    stop(
      "order needs to contain either orders for all dimensions of x or just orders for the selected margin."
    )

  # set margins not to be permuted to identity and copy the rest
  o <- seriate(x, method = "identity")
  if (length(order) <  ndim(x)) ### we only have order for specified margins
    for(i in seq(length(order)))
      o[[margin[i]]] <- order[[i]]
  else
    for (i in margin)
      o[[i]] <- order[[i]]

  # expand identity manual permutations (if any)
  for (i in which(sapply(o, .is_identity_permutation)))
    o[[i]] <- ser_permutation_vector(seq(dim(x)[i]))

  # check
  .check_matrix_perm(x, o)

  perm <- lapply(o, get_order)
  do.call("[", c(list(x), perm, drop = FALSE))
}

.permute_1d <- function(x, order, ...) {
  .nodots(...)

  if (!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  if (length(order) != 1)
    stop("dimensions do not match!")

  if (.is_identity_permutation(order[[1]]))
    return(x)

  perm <- get_order(order, 1)
  if (length(x) != length(perm))
    stop("some permutation vectors do not fit dimension of data!")

  x[perm]
}


# if we used proxy we would say:
#.rearrange_dist <- function (x, order) x[[order]]

.rearrange_dist <- function (x, order) {
  # make C call
  mode(x) <- "double"
  # as.dist seems to make Size numeric and not integer!
  attr(x, "Size") <- as.integer(attr(x, "Size"))
  mode(order) <- "integer"

  d <- .Call("reorder_dist", x, order)

  labels <- if (is.null(labels(x)))
    NULL
  else
    labels(x)[order]

  structure(
    d,
    class   = "dist",
    Size    = length(order),
    Labels  = labels,
    Diag    = FALSE,
    Upper   = FALSE,
    method  = attr(x, "method")
  )
}
