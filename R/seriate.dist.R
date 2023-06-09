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

#' @rdname seriate
#' @export
seriate.dist <-
  function(x,
    method = "Spectral",
    control = NULL,
    ...) {

    ## add ... to control
    control <- c(control, list(...))

    ## check x
    if (anyNA(x)) stop("NAs not allowed in distance matrix x!")
    if (any(x < 0)) stop("Negative distances not supported!")

    if (!is.character(method) || (length(method) != 1L))
      stop("Argument 'method' must be a character string.")
    method <- get_seriation_method("dist", method)

    if (!is.null(control$verbose) &&
        control$verbose)
      cat("Using seriation method: ", method$name, "\n",
        method$description, "\n\n", sep = "")

    tm <- system.time(order <- method$fun(x, control = control))
    if (is.integer(order)) names(order) <- labels(x)[order]

    if (!is.null(control$verbose) &&
        control$verbose)
      cat("Seriation took", tm[1] + tm[2], "sec\n\n")

    ser_permutation(ser_permutation_vector(order, method = method$name))
  }



