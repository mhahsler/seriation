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

.seriate_mean_control <- list(
  transformation = NULL
  )
attr(.seriate_mean_control, "help") <- list(
  transformation = "transformation function applied before calculating means (e.g., scale)"
)

seriate_matrix_mean <- function(x, control = NULL, margin = NULL) {
  control <- .get_parameters(control, .seriate_mean_control)

  if(!is.null(control$transformation))
    x <- control$transformation(x)

  if (1L %in% margin)
    row <- order(rowMeans(x, na.rm = TRUE))
   else
    row <- NA

  if (2L %in% margin)
    col <- order(colMeans(x, na.rm = TRUE))
   else
    col <- NA

  list(row = row, col = col)
}


set_seriation_method(
  "matrix",
  "Mean",
  seriate_matrix_mean,
  "Reorders rows and columns by row and column means.",
  .seriate_mean_control
)
