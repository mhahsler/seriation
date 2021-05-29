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

## image method that makes a proper image plot of a matrix.
## the rows and columns are swapped and the order of the
## columns (original rows) is reversed.


ggpimage <-
  function(x,
    order = NULL,
    axes = "auto",
    upper.tri = TRUE,
    lower.tri = TRUE)
    UseMethod("ggpimage")

### Note for matrix large values are dark, for dist large values are light!
ggpimage.matrix <-
  function(x,
    order = NULL,
    axes = "auto",
    upper.tri = TRUE,
    lower.tri = TRUE) {
    check_installed("ggplot2")

    # check data
    x <- as.matrix(x)
    if (all(is.na(x)))
      stop("all data missing in x.")
    if (any(is.infinite(x)))
      stop("x contains infinite entries.")

    # reorder
    if (!is.null(order))
      x <- permute(x, order)

    # mask triangles
    if (any(!upper.tri ||
        !lower.tri) &&
        nrow(x) != ncol(x))
      stop("Upper or lower triangle can only be suppressed for square matrices!")
    if (!upper.tri)
      x[upper.tri(x)] <- NA
    if (!lower.tri)
      x[lower.tri(x)] <- NA

    # convert to data.frame
    x_df <- data.frame(
      row = factor(rep(seq(nrow(
        x
      )), times = ncol(x)), levels = seq(nrow(x), 1)),
      col = factor(rep(seq(ncol(
        x
      )), each = nrow(x)), levels = seq(ncol(x))),
      x = as.vector(x)
    )

    if (!is.null(rownames(x)))
      levels(x_df[["row"]]) <- rownames(x)
    if (!is.null(colnames(x)))
      levels(x_df[["col"]]) <- colnames(x)

    m <-
      match.arg(axes, c("auto", "x", "y", "rows", "cols", "both", "none"))

    row_labels <- col_labels <- FALSE

    if (m == "auto") {
      row_labels <- nrow(x) >= 25
      col_labels <- ncol(x) >= 25
    }
    else if (m == "x" || m == "cols")
      col_labels <- TRUE
    else if (m == "y" || m == "rows")
      row_labels <- TRUE
    else if (m == "both")
      row_labels <- col_labels <- TRUE

    g <- ggplot2::ggplot(x_df,
      ggplot2::aes(y = row,
        x = col,
        fill = x)) +
      ggplot2::geom_raster() +
      ggplot2::scale_x_discrete(breaks = if (col_labels) ggplot2::waiver() else NULL,
        expand = c(0, 0)) +
      ggplot2::scale_y_discrete(breaks = if (row_labels) ggplot2::waiver() else NULL,
        expand = c(0, 0)) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), axis.title = ggplot2::element_blank())

    # colors for logical
    if (is.logical(x))
      g <- g + ggplot2::scale_fill_manual(values = c("gray", "black"))

    # colors for diverging
    if (any(x < 0, na.rm = TRUE) && any(x > 0, na.rm = TRUE))
      g <- g + ggplot2::scale_fill_gradient2(low = "blue" , mid = "white", high = "red", midpoint = 0)

    g

  }

ggpimage.default <- ggpimage.matrix

## small values are dark
ggpimage.dist <-
  function(x,
    order = NULL,
    axes = "auto",
    upper.tri = TRUE,
    lower.tri = TRUE) {
    check_installed("ggplot2")

    # reorder
    if (!is.null(order))
      x <- permute(x, order)
    order <- NULL

    ggpimage.matrix(as.matrix(x), order, axes, upper.tri, lower.tri) +
      ggplot2::scale_fill_gradient(low = "black",
        high = "white",
        na.value = "white") +
      ggplot2::theme(aspect.ratio = 1)
  }
