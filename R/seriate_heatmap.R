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
  dist_fun = list(row = dist, col = dist),
  seriation_method = list(row = "OLO", col = "OLO"),
  seriation_control = list(row = NULL, col = NULL),
  scale = "none",
  verbose = FALSE
)

seriate_matrix_heatmap <- function(x, control = NULL) {
  control <- .get_parameters(control, .heatmap_contr)

  if (length(control$dist_fun) == 1L)
    control$dist_fun <-
      list(row = control$dist_fun,
        col = control$dist_fun)
  if (length(control$seriation_method) == 1L)
    control$seriation_method <-
      list(row = control$seriation_method,
        col = control$seriation_method)
  if (length(control$seriation_control) == 1L)
    control$seriation_control <-
      list(row = control$seriation_control,
        col = control$seriation_control)

  if (!is.null(control$scale)) {
    if (control$scale == "row")
      x <- t(scale(t(x)))
    if (control$scale == "col")
      x <- scale(x)
  }

  d <- control$dist_fun$row(x)
  o_row <- seriate(d,
    method = control$seriation_method$row,
    control = control$seriation_control$row)

  d <- control$dist_fun$col(t(x))
  o_col <- seriate(d,
    method = control$seriation_method$col,
    control = control$seriation_control$col)

  #names(row) <- rownames(x)[get_order(o_row)]
  #names(col) <- colnames(x)[get_order(o_col)]

  list(row = o_row, col = o_col)
}

set_seriation_method(
  "matrix",
  "Heatmap",
  seriate_matrix_heatmap,
  "Calculate distances for row and column vectors, and seriate. If only a single distance function or seriation method is specified, then it is used for rows and columns. The default seriation method is optimal leaf ordering (OLO) which  perform hierarchical clustering and reorder the dentrograms.",
  .heatmap_contr
)
