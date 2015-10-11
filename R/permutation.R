#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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


######################################################
## permutations

## constructor
ser_permutation <- function(x, ...) {
  x <- c(list(x), list(...))
  
  x <- lapply(x, FUN = function(obj) {
    if(is(obj, "ser_permutation")) return(obj)
    if(is(obj, "ser_permutation_vector")) return(list(obj))
    return(list(ser_permutation_vector(obj)))
    })
  
  x <- unlist(x, recursive = FALSE)
  class(x) <- c("ser_permutation", "list")
  x
}

## so we can say get_order to permutations
get_order.ser_permutation <- function(x, dim = 1, ...) get_order(x[[dim]])

## print et al
print.ser_permutation <- function(x, ...) {
  writeLines(c(
    gettextf("object of class %s",
      paste(sQuote(class(x)), collapse = ", ")),
    gettextf("contains permutation vectors for %d-mode data\n",
      length(x))
  ))
  
  print(data.frame("vector length" = sapply(x, 
    FUN = function(o) if(.is_identity_permutation(o)) NA_integer_ else length(o)),
    "seriation method" =
      sapply(x, get_method, printable = TRUE),
    check.names = FALSE))
  
  invisible(x)
}

## fake summary (we dont really provide a summary, 
## but summary produces now a reasonable result --- same as print)
summary.ser_permutation <- function(object, ...) 
  object

c.ser_permutation <- function(..., recursive = FALSE) 
  do.call("ser_permutation", list(...)) 

## fixme [[<- needs to check for ser_permutation_vector

"[.ser_permutation" <- function(object, i, ...) 
  do.call("ser_permutation", unclass(object)[i])
