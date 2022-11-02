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

## ser_permutation_vector represents a single permutation represented as an
## integer vector or a hclust object.

## Constructor
## x can be
##  * an integer vector
##  * a hclust or dendrogram object (leaf order)
##  * NA represents the identity permutation
##  * a ser_permutation (list) of length 1

#' Class ser_permutation_vector -- A Single Permutation Vector for Seriation
#'
#' The class `ser_permutation_vector`
#' represents a single permutation vector.
#'
#' A permutation vector
#' maps a set of \eqn{n} objects \eqn{\{O_1, O_2, ..., O_n\}}{{O_1, O_2, ..., O_n}} onto itself.
#'
#' __Ordering Representation:__
#' In \pkg{seriation} we represent a permutation \eqn{\pi}{\pi}
#' as a vector which lists the objects' indices in their permuted order. This can
#' be seen as replacing the object in position \eqn{i} with the object
#' in position \eqn{\pi(i)}.
#' For example, the permutation vector \eqn{\langle3, 1, 2\rangle}{<3, 1, 2>} indicates that in
#' first position is the object with index 3 then the object with index 1 and finally
#' the object with index 2. This representation is often called a (re)arrangement or ordering.
#' The ordering can be extracted from a permutation vector object
#' via [get_order()]. Such an ordering can be directly used
#' to subset the list of original objects with `"["` to apply the permutation.
#'
#' __Rank Representation:__
#' An alternative way to specify a permutation is via a list of the ranks
#' of the objects after permutation. This representation is often called
#' a map or substitution. Ranks can be extracted from a permutation vector using [get_rank()].
#'
#' __Permutation Matrix:__
#' Another popular representation is a permutation matrix which performs
#' permutations using matrix multiplication. A permutation matrix can be obtained
#' using [get_permutation_matrix()].
#'
#' `ser_permutation_vector` objects are usually packed into
#' a [ser_permutation] object
#' which is a collection (a `list`) of \eqn{k} permutation vectors for \eqn{k}-mode data.
#'
#' The constructor `ser_permutation_vector()`
#' checks if the permutation vector is valid
#' (i.e. if all integers occur exactly once).
#'
#' @family permutation
#'
#' @param x,object an object which contains a permutation vector (currently an
#'     integer vector or an object of class [hclust]). The value `NA`
#'     creates an identity permutation.
#' @param method a string representing the method used to obtain the
#'     permutation vector.
#' @param ... further arguments.
#'
#' @returns  An object of class `ser_permutation_vector`.
#' @author Michael Hahsler
#'
#' @examples
#' p <- ser_permutation_vector(sample(10), "random")
#' p
#'
#' ## some methods
#' length(p)
#' get_method(p)
#' get_order(p)
#' get_rank(p)
#' get_permutation_matrix(p)
#'
#' r <- rev(p)
#' r
#' get_order(r)
#'
#' ## create a identity permutation vector (with unknown length)
#' ip <- ser_permutation_vector(NA)
#' ip
#'
#' @keywords classes
#' @export
ser_permutation_vector <- function(x, method = NULL) {
  if (inherits(x, "ser_permutation_vector"))
    return(x)

  if (inherits(x, "hclust") || inherits(x, "dendrogram")) {
    # nothing to do
  } else if (length(x) == 1 && is.na(x)) {
    x <- NA_integer_
    attr(x, "method") <- "identity permutation"
  } else if (is.numeric(x)) {
    nm <- names(x)
    x <- as.integer(x)
    names(x) <- nm
  } else if (inherits(x, "ser_permutation") && length(x) == 1) {
    x <- x[[1]]
  } else {
    stop("x does not contain a supported permutation.")
  }

  if (!is.null(method))
    attr(x, "method") <- method

  class(x) <- c("ser_permutation_vector", class(x))
  .valid_permutation_vector(x)
  x
}

#' @rdname ser_permutation_vector
#' @param recursive ignored
#' @export
c.ser_permutation_vector <- function(..., recursive = FALSE)
  do.call("ser_permutation", list(...))




## reverse
#' @rdname ser_permutation_vector
#' @export
rev.ser_permutation_vector <- function(x) {
  if (inherits(x, "hclust")) {
    x$order <- rev(x$order)
    x
  }
  else
    ser_permutation_vector(rev(get_order(x)), method = get_method(x))
}


## currently method is an attribute of permutation
#' @rdname ser_permutation_vector
#' @param printable a logical; prints "unknown" instead of `NULL` for non-existing methods.
#' @export
get_method <- function(x, printable = FALSE) {
  method <- attr(x, "method")

  if (printable && is.null(method))
    method <- "unknown"
  method
}


## print et al
#' @rdname ser_permutation_vector
#' @export
length.ser_permutation_vector <- function(x) {
  if (!.is_identity_permutation(x))
    length(get_order(x))
  else
    0L
}

#' @rdname ser_permutation_vector
#' @export
print.ser_permutation_vector <-
  function(x, ...)
  {
    writeLines(c(
      gettextf("object of class %s",
        paste(sQuote(class(
          x
        )), collapse = ", ")),
      gettextf("contains a permutation vector of length %d", length(x)),
      gettextf("used seriation method: '%s'",
        get_method(x, printable = TRUE))
    ))
    invisible(x)
  }

## fake summary (we don't really provide a summary,
## but summary produces now a reasonable result --- same as print)
#' @rdname ser_permutation_vector
#' @export
summary.ser_permutation_vector <- function(object, ...) {
  object
}


## helpers

## an identity permutation is a single NA.
.is_identity_permutation <- function(x) is.na(x[1])

## calls stop if the vector is not valid
.valid_permutation_vector <- function(x) {

  ## identity vector is always valid
  if (.is_identity_permutation(x))
    return(invisible(TRUE))

  ## valid permutations have a get_order function implemented
  perm <- get_order(x)
  valid <- TRUE

  tab <- table(perm)
  if (any(tab != 1))
    valid <- FALSE
  if (length(tab) != length(perm)
    || any(names(tab) != sequence(length(perm))))
    valid <- FALSE

  if (!valid)
    stop("Invalid permutation vector!\nVector: ",
      paste(perm, collapse = ", "))

  invisible(valid)
}

.valid_permutation_matrix <- function(x) {
  if (any(rowSums(x) != 1) || any(colSums(x) != 1) ||
      any(x != 1 & x != 0))
    stop("Not a valid permutation matrix")

  invisible(TRUE)
}
