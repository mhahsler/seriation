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

#' @rdname pimage
#' @include pimage.R
#' @export
ggpimage <- function(x,
  order = NULL,
  zlim = NULL,
  upper_tri = TRUE,
  lower_tri = TRUE,
  diag = TRUE,
  row_labels = NULL,
  col_labels = NULL,
  prop = TRUE,
  flip_axes = FALSE,
  reverse_columns = FALSE)
  UseMethod("ggpimage")


### Note for matrix large values are dark, for dist large values are light!
#' @rdname pimage
#' @export
ggpimage.matrix <- function(x,
  order = NULL,
  zlim = NULL,
  upper_tri = TRUE,
  lower_tri = TRUE,
  diag = TRUE,
  row_labels = NULL,
  col_labels = NULL,
  prop = TRUE,
  flip_axes = FALSE,
  reverse_columns = FALSE)
{
  check_installed("ggplot2")

  x <- as.matrix(x)

  # check data
  if (all(is.na(x)))
    stop("all data missing in x.")
  if (any(is.infinite(x)))
    stop("x contains infinite entries.")

  # reorder
  if (!is.null(order))
    x <- permute(x, order)

  # mask triangles
  if (any(!upper_tri ||
      !lower_tri ||
      !diag) &&
      nrow(x) != ncol(x))
    stop("Upper triangle, lower triangle or diag can only be suppressed for square matrices!")
  if (!upper_tri)
    x[upper.tri(x)] <- NA
  if (!lower_tri)
    x[lower.tri(x)] <- NA
  if (!diag)
    diag(x) <- NA


  # reverse order of columns
  if (reverse_columns)
    x <- x[, seq(ncol(x), 1)]

  # change x and y?
  if (flip_axes) {
    x <- t(x)
    tmp <- row_labels
    row_labels <- col_labels
    col_labels <- tmp
  }


  # plot
  g <-
    .ggpimage_empty(
      x,
      zlim = zlim,
      row_labels = row_labels,
      col_labels = col_labels,
      prop = prop,
      expand = FALSE
    )

  g <- g + ggplot2::geom_raster(ggplot2::aes(fill = x))
  g
}

ggpimage.default <- ggpimage.matrix

## small values are dark
#' @rdname pimage
#' @export
ggpimage.dist <-
  function(x,
    order = NULL,
    zlim = NULL,
    upper_tri = FALSE,
    lower_tri = TRUE,
    diag = FALSE,
    row_labels = NULL,
    col_labels = NULL,
    prop = TRUE,
    flip_axes = FALSE,
    reverse_columns = FALSE) {
    check_installed("ggplot2")

    # reorder specific for dist (we have only a single permutation)
    if (!is.null(order))
      x <- permute(x, order)

    if (flip_axes)
      warning("flipping axes has no effect for distance matrices.")

    g <- ggpimage.matrix(
      as.matrix(x),
      order = NULL,
      zlim = zlim,
      upper_tri,
      lower_tri,
      diag,
      row_labels,
      col_labels,
      prop = prop,
      flip_axes = FALSE,
      reverse_columns = reverse_columns
    )

    # reverse color for dist
    suppressMessages(g <-
        g + .gg_sequential_pal(dist = TRUE, limits = zlim))

    g
  }

### Note for matrix large values are dark, for dist large values are light!
.ggpimage_empty <- function(x,
  zlim = NULL,
  row_labels = NULL,
  col_labels = NULL,
  prop = TRUE,
  expand = TRUE) {
  check_installed("ggplot2")

  x <- as.matrix(x)

  # check data
  if (all(is.na(x)))
    stop("all data missing in x.")
  if (any(is.infinite(x)))
    stop("x contains infinite entries.")

  # deal with row/col labels
  if (!is.null(row_labels) && !is.logical(row_labels)) {
    if (length(row_labels) != nrow(x))
      stop("Length of row_labels does not match the number of rows of x.")
    rownames(x) <- row_labels
    row_labels <- TRUE
  }

  if (!is.null(col_labels) && !is.logical(col_labels)) {
    if (length(col_labels) != ncol(x))
      stop("Length of col_labels does not match the number of columns of x.")
    colnames(x) <- col_labels
    col_labels <- TRUE
  }

  if (is.null(row_labels))
    if (!is.null(rownames(x)) &&
        nrow(x) < 25) {
      row_labels <- TRUE
    } else{
      row_labels <- FALSE
    }

  if (is.null(col_labels))
    if (!is.null(colnames(x)) &&
        ncol(x) < 25) {
      col_labels <- TRUE
    } else{
      col_labels <- FALSE
    }

  if (is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if (is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))

  # convert to data.frame with row, col and x
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
    levels(x_df[["row"]]) <- rev(rownames(x))
  if (!is.null(colnames(x)))
    levels(x_df[["col"]]) <- colnames(x)

  # plot
  g <- ggplot2::ggplot(x_df,
    ggplot2::aes(y = row, x = col))

  # axes (row and col labels)
  if (expand)
    expand <- ggplot2::waiver()
  else
    expand <- c(0, 0)

  if (col_labels)
    breaksCol <- ggplot2::waiver()
  else
    breaksCol <- NULL
  if (row_labels)
    breaksRow <- ggplot2::waiver()
  else
    breaksRow <- NULL

  g <- g +
    ggplot2::scale_x_discrete(breaks = breaksCol, expand = expand) +
    ggplot2::scale_y_discrete(breaks = breaksRow, expand = expand)

  # no axis or legend labels
  g <- g +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL)

  g <-
    g + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    ))

  if (prop)
    g <- g + ggplot2::theme(aspect.ratio = nrow(x) / ncol(x))

  # colors scalesi
  if (is.logical(x)) {
    g <-
      g + .gg_logical_pal()

    # colors for diverging
  } else if (any(x < 0, na.rm = TRUE) && any(x > 0, na.rm = TRUE)) {
    g <-
      g + .gg_diverge_pal(limits = zlim)

  } else {
    g <-
      g + .gg_sequential_pal(limits = zlim)
  }

  g
}


