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
  seriation_method = list(row = "OLO_complete", col = "OLO_complete"),
  seriation_control = list(row = NULL, col = NULL),
  scale = "none",
  verbose = FALSE
)

 attr(.heatmap_contr, "help") <- list(
  dist_fun = "A named list with functions to calulate row and column distances",
  seriation_method = "A named list with row and column seriation methods",
  seriation_control = "named list with control parameters for the seriation methods",
  scale = 'Scale "rows", "cols", or "none"'
)


seriate_matrix_heatmap <-
  function(x,
           control = NULL,
           margin = seq_along(dim(x))) {
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
      if (control$scale == "rows")
        x <- t(scale(t(x)))
      if (control$scale == "cols")
        x <- scale(x)
    }

    if (1L %in% margin) {
      d <- control$dist_fun$row(x)

      if (tolower(control$seriation_method$row) == "hc_mean")
        o_row <- ser_permutation_vector(seriate_hc_mean(d, x, control$seriation_control$row), method = "HC_Mean")
      else
        o_row <- seriate(
          d,
          method = control$seriation_method$row,
          control = control$seriation_control$row
        )[[1]]
    } else
      o_row <- NA

    if (2L %in% margin) {
      x <- t(x)
      d <- control$dist_fun$col(x)

      if (tolower(control$seriation_method$col) == "hc_mean")
        o_col <- ser_permutation_vector(seriate_hc_mean(d, x, control$seriation_control$col), method = "HC_Mean")
      else
        o_col <- seriate(
          d,
          method = control$seriation_method$col,
          control = control$seriation_control$col
        )[[1]]
    } else
      o_col <- NA

    #names(row) <- rownames(x)[get_order(o_row)]
    #names(col) <- colnames(x)[get_order(o_col)]

    list(row = o_row, col = o_col)
  }


seriate_hc_mean <- function(d, x, control = NULL) {
  if (missing(x))
    stop("data matrix x needs to be specified for leaf order with mean reordering.")

  hc <- stats::as.hclust(stats::reorder(
    stats::as.dendrogram(seriate_dist_hc(d, control)),
    wts = rowSums(x, na.rm = TRUE)
  ))
  hc$call <- match.call()
  hc$method <- "hclust + mean reordering"
  hc$dist.method <- attr(d, "method")

  hc
}

set_seriation_method(
  "matrix",
  "Heatmap",
  seriate_matrix_heatmap,
  "Calculates distances for rows and columns and then independently applies the specified seriation method for distances.",
  control = .heatmap_contr
)
