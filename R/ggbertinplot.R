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


ggbertinplot <- function(x,
  order = NULL,
  highlight = TRUE,
  labRow = TRUE,
  labCol = TRUE,
  reverse = FALSE,
  prop = FALSE,
  geom = "bar") {
  check_installed("ggplot2")

  if (!is.matrix(x))
    stop("Argument 'x' must be a matrix.")

  if (is.logical(highlight) && highlight)
    highlight <- mean(x, na.rm = TRUE)

  # TODO: rename to circle, box etc.
  # Add square
  geom <-
    match.arg(tolower(geom),
      choices = c("raster", "block", "circle", "line", "bar", "none"))

  # reorder
  if (!is.null(order))
    x <- permute(x, order)

  # change x and y?
  if (!reverse) {
    x <- t(x)
    tmp <- labRow
    labRow <- labCol
    labCol <- tmp
  }

  g <-
    .ggpimage_empty(
      x,
      labRow = labRow,
      labCol = labCol,
      prop = prop,
      expand = geom != "raster"
    )

  # put col labels on top (message about replacing scale for x)
  suppressMessages(
    g <- g +
      ggplot2::scale_x_discrete(position = "top", expand = if (geom != "raster") ggplot2::waiver() else c(0, 0)) +
      ggplot2::scale_y_discrete(position = "right", expand = if (geom != "raster") ggplot2::waiver() else c(0, 0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0, vjust = .5)) +
      ggplot2::theme(legend.position="bottom")
  )

  # add geom

  # raster does not use highlight
  if (geom == "raster")
    g <- g + ggplot2::geom_raster(ggplot2::aes(fill = x))

  if (geom == "circle")
    if (highlight) {
      g <-
        g + ggplot2::geom_point(ggplot2::aes(size = x,
          fill = x > highlight),
          color = "black",
          pch = 21) + ggplot2::scale_fill_manual(values = c("white", "black")) +
        ggplot2::guides(fill = FALSE, size = FALSE)
    } else{
      g <- g + ggplot2::geom_point(ggplot2::aes(size = x))
    }

  if (geom == "block")
    if (highlight) {
      g <-
        g + ggplot2::geom_tile(
          ggplot2::aes(
            x = col,
            y = row,
            height = x / max(x) * .8,
            fill = x > highlight
          ),
          width = .8,
          color = "black"
        ) +
        ggplot2::scale_fill_manual(values = c("white", "black")) +
        ggplot2::guides(fill = FALSE)
    } else{
      g <- g +
        ggplot2::geom_tile(ggplot2::aes(height = x / max(x) * .9), width = .8)
    }

  # no highlight for line
  if (geom == "line")
    g <- g +
    ggplot2::geom_line(ggplot2::aes(x = col, y = x, group = row)) +
    # Note: facets display the lowest level first so we need to reverse them
    ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(row, rev(as.integer(row))))) +
    ggplot2::theme(
      strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
      strip.background = ggplot2::element_blank()
    )

  if (geom == "bar")
    if (highlight) {
      g <- g +
        ggplot2::geom_bar(ggplot2::aes(
          x = col,
          y = x,
          group = row,
          fill = x > highlight
        ),
          stat = "identity", color = "black", width = .8) +
        # Note: facets display the lowest level first so we need to reverse them
        ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(row, rev(as.integer(row))))) +
        ggplot2::theme(
          strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
          strip.background = ggplot2::element_blank()
        ) +
        ggplot2::scale_fill_manual(values = c("white", "black")) +
        ggplot2::guides(fill = FALSE)
    } else{
      g <- g +
        ggplot2::geom_bar(ggplot2::aes(x = col,
          y = x,
          group = row),
          stat = "identity", width = .8) +
        # Note: facets display the lowest level first so we need to reverse them
        ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(row, rev(as.integer(row))))) +
        ggplot2::theme(
          strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
          strip.background = ggplot2::element_blank()
        )
    }

  g
}
