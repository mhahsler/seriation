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

#' Class ser_permutation -- A Collection of Permutation Vectors for Seriation
#'
#' The class `ser_permutation` is a collection of permutation vectors
#' (see class [ser_permutation_vector]), one for each dimension (mode)
#' of the data to be permuted.
#'
#' @family permutation
#'
#' @param x,object an object of class `ser_permutation_vector` or
#'     any object which can be converted into
#'     a object of class `ser_permutation` (e.g. an integer
#'       vector).
#' @param ... vectors for further dimensions.
#'
#' @returns An object of class `ser_permutation`.
#'
#' @author Michael Hahsler
#' @examples
#' o <- ser_permutation(1:5, 10:1)
#' o
#'
#' ## length (number of dimensions)
#' length(o)
#'
#' ## get permutation vector for 2nd dimension
#' get_order(o, 2)
#'
#' ## reverse dimensions
#' o[2:1]
#'
#' ## combine
#' o <- c(o, ser_permutation(1:15))
#' o
#'
#' ## get an individual permutation
#' o[[2]]
#'
#' ## reverse the order of a permutation
#' o[[2]] <- rev(o[[2]])
#' get_order(o,2)
#' @keywords classes
#' @export
ser_permutation <- function(x, ...) {
  x <- c(list(x), list(...))

  x <- lapply(
    x,
    FUN = function(obj) {
      if (inherits(obj, "ser_permutation"))
        return(obj)
      if (inherits(obj, "ser_permutation_vector"))
        return(list(obj))
      return(list(ser_permutation_vector(obj)))
    }
  )

  x <- unlist(x, recursive = FALSE)
  class(x) <- c("ser_permutation", "list")
  x
}

#' @rdname ser_permutation
#' @export
print.ser_permutation <- function(x, ...) {
  writeLines(c(
    gettextf("object of class %s",
      paste(sQuote(class(
        x
      )), collapse = ", ")),
    gettextf("contains permutation vectors for %d-mode data\n",
      length(x))
  ))

  print(
    data.frame(
      "vector length" = sapply(
        x,
        FUN = function(o)
          if (.is_identity_permutation(o))
            NA_integer_
        else
          length(o)
      ),
      "seriation method" =
        sapply(x, get_method, printable = TRUE),
      check.names = FALSE
    )
  )

  invisible(x)
}

## fake summary (we don't really provide a summary,
## but summary produces now a reasonable result --- same as print)
#' @rdname ser_permutation
#' @export
summary.ser_permutation <- function(object, ...)
  object

#' @rdname ser_permutation
#' @param recursive ignored.
#' @export
c.ser_permutation <- function(..., recursive = FALSE)
  do.call("ser_permutation", list(...))

## fixme [[<- needs to check for ser_permutation_vector

#' @rdname ser_permutation
#' @param i index of the dimension(s) to extract.
#' @export
"[.ser_permutation" <- function(object, i, ...)
  do.call("ser_permutation", unclass(object)[i])

is.ser_permutation <- function(x)
  inherits(x, "ser_permutation") | inherits(x, "ser_permutation_vector")
