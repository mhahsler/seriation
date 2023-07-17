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

.lle_contr <- list(
  k = 30,
  reg = 2
)

attr(.lle_contr, "help") <- list(
  k = "used number of neighbors",
  reg = "regularization method (see ? lle)"
)

seriate_lle <- function(x, control = NULL, margin) {
  param <- .get_parameters(control, .lle_contr)

  o <- list(row = NA, col = NA)

  if (1L %in% margin) {
    score <- lle(x, m = 1, k = param$k, reg = param$reg)
    os <- order(score)
    o$row <- structure(os, names = rownames(x)[os], configuration = score)
  }

  if (2L %in% margin) {
    x <- t(x)
    score <- lle(x, m = 1, k = param$k, reg = param$reg)
    os <- order(score)
    o$col <- structure(os, names = rownames(x)[os], configuration = score)
  }

  o
}

set_seriation_method(
  "matrix",
  "LLE",
  seriate_lle,
  "Find an order using 1D locally linear embedding.\n",
  .lle_contr,
  randomized = FALSE
)
