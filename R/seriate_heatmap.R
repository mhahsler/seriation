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


## calculate distances for rows and columns, perform hclust and reorder.
.heatmap_contr <- list(
  distfun = dist,
  method = "OLO",
  control = NULL,
  scale = c("none"),
  verbose = FALSE
)

seriate_matrix_heatmap <- function(x, control = NULL) {
  control <- .get_parameters(control, .heatmap_contr)

  control$scale <- match.arg(control$scale, choices = c("none", "row", "column"))

  if (!is.null(control$scale)) {
    if (control$scale == "row")
      x <- t(scale(t(x)))
    if (control$scale == "col")
      x <- scale(x)
  }

  dist_row <- control$distfun(x)
  o_row <- seriate(dist_row,
    method = control$method,
    control = control$control)[[1]]

  dist_col <- control$distfun(t(x))
  o_col <- seriate(dist_col,
    method = control$method,
    control = control$control)[[1]]

  #names(row) <- rownames(x)[get_order(o_row)]
  #names(col) <- colnames(x)[get_order(o_col)]

  list(row = o_row, col = o_col)
}

set_seriation_method(
  "matrix",
  "Heatmap",
  seriate_matrix_heatmap,
  "Calculate distances for row and column vectors, perform hierarchical clustering and reorder the dentrograms.",
  .heatmap_contr
)
