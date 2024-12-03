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

## seriate general arrays

.seriate_array_helper <- function(x,
  method = "PCA",
  control = NULL,
  margin = seq(ndim(x)),
  datatype = "array",
  ...) {
  ## add ... to control
  if (any(!margin %in% seq(ndim(x))))
    stop("illegal margin specified.")

  control <- c(control, list(...))

  if (!is.character(method) || (length(method) != 1L))
    stop("Argument 'method' must be a character string.")

  if (any(dim(x) == 0L))
    stop("All dimensions need to have at least one object.")

  method <- get_seriation_method(datatype, method)

  if (!is.null(control$verbose) &&
      control$verbose)
    cat("Using seriation method: ", method$name, "\n",
        method$description, "\n\n", sep = "")

  tm <- system.time(order <- method$fun(x, control, margin))

  if (!is.null(control$verbose) &&
      control$verbose)
    cat("Seriation took", tm[1] + tm[2], "sec\n\n")

  for (i in margin)
    if (!is.null(dimnames(x)[[i]]) &&
        is.integer(order[[i]]))
      names(order[[i]]) <- dimnames(x)[[i]][order[[i]]]
  perm <- do.call("ser_permutation",
    unname(lapply(
      order, "ser_permutation_vector", method$name
    )))

  ### make non-seriated margins identity permutations
  rem <- which(!seq(ndim(x)) %in% margin)
  if (length(rem) > 0) {
    perm_ident <- seriate(x, method = "Identity")
    perm[[rem]] <- perm_ident[[rem]]
  }

  perm
}

#' @rdname seriate
#' @include seriate.matrix.R
#' @export
seriate.array <- function(x,
  method = "PCA",
  control = NULL,
  margin = seq(length(dim(x))),
  rep = 1L,
  ...) {
  if (rep > 1L)
    return(seriate_rep(x, method, control, rep = rep, margin = margin, ...))

  .seriate_array_helper(x,
    method,
    control,
    margin,
    datatype = "array",
    ...)
}
