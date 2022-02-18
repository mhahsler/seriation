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

#' @rdname bertinplot
#' @export
ggbertinplot <- function(x,
  order = NULL,
  geom = "bar",
  highlight = TRUE,
  row_labels = TRUE,
  col_labels = TRUE,
  flip_axes = TRUE,
  prop = FALSE,
  ...) {
  check_installed("ggplot2")

  if (!is.matrix(x))
    stop("Argument 'x' must be a matrix.")

  geom <-
    match.arg(tolower(geom),
      choices = c("tile", "rectangle", "circle", "line", "bar", "none"))

  # reorder
  if (!is.null(order))
    x <- permute(x, order)

  # change x and y?
  if (flip_axes) {
    x <- t(x)
    tmp <- row_labels
    row_labels <- col_labels
    col_labels <- tmp
  }

  if (is.logical(highlight) && highlight)
    highlight <- mean(x, na.rm = TRUE)

  g <-
    .ggpimage_empty(
      x,
      row_labels = row_labels,
      col_labels = col_labels,
      prop = prop,
      expand = geom != "raster"
    )

  if (col_labels)
    breaksCol <- ggplot2::waiver()
  else
    breaksCol <- NULL
  if (row_labels)
    breaksRow <- ggplot2::waiver()
  else
    breaksRow <- NULL

  # put col labels on top (message about replacing scale for x)
  suppressMessages(
    g <- g +
      ggplot2::scale_x_discrete(
        breaks = breaksRow,
        position = "top",
        expand = if (geom != "raster")
          ggplot2::waiver()
        else
          c(0, 0)
      ) +
      ggplot2::scale_y_discrete(
        breaks = breaksCol,
        position = "right",
        expand = if (geom != "raster")
          ggplot2::waiver()
        else
          c(0, 0)
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 0, vjust = .5)) +
      ggplot2::theme(legend.position = "bottom")
  )

  # add geom

  # raster does not use highlight
  if (geom == "tile")
    g <- g + ggplot2::geom_raster(ggplot2::aes(fill = x))

  if (geom == "circle")
    if (highlight) {
      suppressMessages(
        g <-
          g + ggplot2::geom_point(
            ggplot2::aes(size = x,
              fill = x > highlight),
            color = "black",
            pch = 21
          ) + .gg_logical_pal() +
          ggplot2::guides(fill = "none", size = "none")
      )
    } else{
      g <- g + ggplot2::geom_point(ggplot2::aes(size = x))
    }

  if (geom == "rectangle")
    if (highlight) {
      suppressMessages(
        g <-
          g + ggplot2::geom_tile(
            ggplot2::aes(
              x = col,
              y = row,
              height = x / max(x, na.rm = TRUE) * .8,
              width = x / max(x, na.rm = TRUE) * .8,
              fill = x > highlight
            ),
            color = "black"
          ) +
          .gg_logical_pal() +
          ggplot2::guides(fill = "none")
      )
    } else{
      g <- g +
        ggplot2::geom_tile(ggplot2::aes(height = x / max(x) * .9), width = .8)
    }

  # TODO: do not display facet labels when row_labels == FALSE

  # no highlight for line
  if (geom == "line")
    g <- g +
    ggplot2::geom_line(ggplot2::aes(x = col, y = x, group = row)) +
    # Note: facets display the lowest level first so we need to reverse them
    ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(row, rev(as.integer(
      row
    ))))) +
    ggplot2::theme(
      strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
      strip.background = ggplot2::element_blank()
    )

  if (geom == "bar")
    if (highlight) {
      suppressMessages(
        g <- g +
          ggplot2::geom_bar(
            ggplot2::aes(
              x = col,
              y = x,
              group = row,
              fill = x > highlight
            ),
            stat = "identity",
            color = "black",
            width = .8
          ) +
          # Note: facets display the lowest level first so we need to reverse them
          ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(
            row, rev(as.integer(row))
          ))) +
          ggplot2::theme(
            strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
            strip.background = ggplot2::element_blank()
          ) +
          .gg_logical_pal() +
          ggplot2::guides(fill = "none")
      )
    } else{
      g <- g +
        ggplot2::geom_bar(ggplot2::aes(x = col,
          y = x,
          group = row),
          stat = "identity",
          width = .8) +
        # Note: facets display the lowest level first so we need to reverse them
        ggplot2::facet_grid(rows = ggplot2::vars(stats::reorder(row, rev(
          as.integer(row)
        )))) +
        ggplot2::theme(
          strip.text.y.right = ggplot2::element_text(angle = 0, color = "black"),
          strip.background = ggplot2::element_blank()
        )
    }

  g
}
