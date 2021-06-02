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

ggpimage <- function(x,
  order = NULL,
  upper.tri = TRUE,
  lower.tri = TRUE,
  labRow = NULL,
  labCol = NULL,
  prop = TRUE,
  flip = FALSE,
  geom = "raster")
  UseMethod("ggpimage")

### Note for matrix large values are dark, for dist large values are light!
ggpimage.matrix <- function(x,
  order = NULL,
  upper.tri = TRUE,
  lower.tri = TRUE,
  labRow = NULL,
  labCol = NULL,
  prop = TRUE,
  flip = FALSE,
  geom = "raster") {
  check_installed("ggplot2")

  x <- as.matrix(x)

  # check data
  if (all(is.na(x)))
    stop("all data missing in x.")
  if (any(is.infinite(x)))
    stop("x contains infinite entries.")

  # check geom
  geom <-
    match.arg(tolower(geom),
      choices = c("raster", "point", "tile", "lines", "none"))

  # change x and y
  if (flip) {
    x <- t(x)
    if (!is.null(order))
      order <- rev(order)
  }

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

  # deal with row/col labels
  if (!is.null(labRow) && !is.logical(labRow)) {
    if (length(labRow) != nrow(x))
      stop("Length of labRow does not match the number of rows of x.")
    rownames(x) <- labRow
    labRow <- TRUE
  }

  if (!is.null(labCol) && !is.logical(labCol)) {
    if (length(labCol) != ncol(x))
      stop("Length of labCol does not match the number of columns of x.")
    colnames(x) <- labCol
    labCol <- TRUE
  }

  if (is.null(labRow))
    if (!is.null(rownames(x)) &&
        nrow(x) < 25) {
      labRow <- TRUE
    } else{
      labRow <- FALSE
    }
  if (is.null(labCol))
    if (!is.null(colnames(x)) &&
        ncol(x) < 25) {
      labCol <- TRUE
    } else{
      labCol <- FALSE
    }

  if (is.null(rownames(x)))
    rownames(x) <- seq(nrow(x))
  if (is.null(colnames(x)))
    colnames(x) <- seq(ncol(x))

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

  # plot
  g <- ggplot2::ggplot(x_df,
    ggplot2::aes(y = row, x = col))

  # axes (row and col labels)
  if (geom == "raster")
    expand <- c(0, 0)
  else
    expand <- ggplot2::waiver()

  if (labCol)
    breaksCol <- ggplot2::waiver()
  else
    breaksCol <- NULL
  if (labRow)
    breaksRow <- ggplot2::waiver()
  else
    breaksRow <- NULL

  g <- g +
    ggplot2::scale_x_discrete(breaks = breaksCol, expand = expand) +
    ggplot2::scale_y_discrete(breaks = breaksRow, expand = expand)

  # no axis or legend labels
  g <- g +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  g <-
    g + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    ))


  if (geom == "raster")
    g <- g + ggplot2::geom_raster(ggplot2::aes(fill = x))

  ### FIXME: circles are clipped!
  if (geom == "point")
    g <- g + ggplot2::geom_point(ggplot2::aes(size = x))

  if (geom == "tile")
    g <- g +
    ggplot2::geom_tile(ggplot2::aes(height = x / max(x) * .9), width = .8)

  if (geom == "lines")
    g <- g +
    ggplot2::geom_line(ggplot2::aes(x = col, y = x, group = row)) +
    ggplot2::facet_grid(rows = ggplot2::vars(row)) +
    ggplot2::theme(
      strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
      strip.background = ggplot2::element_blank(),
    )


  # colors for logical
  if (is.logical(x))
    g <-
    g + ggplot2::scale_fill_manual(values = c("gray", "black"))

  # colors for diverging
  if (any(x < 0, na.rm = TRUE) && any(x > 0, na.rm = TRUE))
    g <-
    g + ggplot2::scale_fill_gradient2(
      low = "blue" ,
      mid = "white",
      high = "red",
      midpoint = 0
    )

  if (prop)
    g <- g + ggplot2::theme(aspect.ratio = nrow(x) / ncol(x))

  g

}

ggpimage.default <- ggpimage.matrix

## small values are dark
ggpimage.dist <-
  function(x,
    order = NULL,
    upper.tri = TRUE,
    lower.tri = TRUE,
    labRow = NULL,
    labCol = NULL,
    prop = TRUE,
    flip = FALSE,
    geom = "raster") {
    check_installed("ggplot2")

    # reorder specific for dist (we have only a single permutation)
    if (!is.null(order))
      x <- permute(x, order)

    if (flip)
      warning("flipping axes has no effect for distance matrices.")

    ggpimage.matrix(
      as.matrix(x),
      order = NULL,
      upper.tri,
      lower.tri,
      labRow,
      labCol,
      prop = prop,
      flip = FALSE,
      geom = geom
    ) +
      ggplot2::scale_fill_gradient(low = "black",
        high = "white",
        na.value = "white")
  }
