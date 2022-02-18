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

#' @rdname dissplot
#' @export
ggdissplot <- function(x,
  labels = NULL,
  method = "Spectral",
  control = NULL,
  lower_tri = TRUE,
  upper_tri = "average",
  diag = TRUE,
  cluster_labels = TRUE,
  cluster_lines = TRUE,
  reverse_columns = FALSE,
  ...) {
  check_installed("ggplot2")

  # make x dist
  if (!inherits(x, "dist")) {
    if (is.matrix(x) && isSymmetric(x))
      x <- as.dist(x)
    else
      stop("Argument 'x' cannot safely be coerced to class 'dist'.")
  }

  x <- .arrange_dissimilarity_matrix(x,
    labels = labels,
    method = method,
    control = control)

  m  <- .average_tri(x,
    lower_tri = lower_tri,
    upper_tri = upper_tri,
    diag = diag)

  k       <- x$k
  dim     <- attr(x$x_reordered, "Size")
  labels  <- x$labels
  labels_unique <- unique(labels)

  # So we can add cluster labels later
  if (cluster_labels)
    colnames(m) <- seq(ncol(m))

  g <- ggpimage(m, reverse_columns = reverse_columns, ...)

  # add cluster lines and labels
  if (!is.null(labels)) {
    cluster_width		<- tabulate(labels)[labels_unique]
    cluster_cuts		<- cumsum(cluster_width)
    cluster_center	<- cluster_cuts - cluster_width / 2

    clusters <-
      data.frame(
        center = cluster_center,
        cut = cluster_cuts,
        width = cluster_width,
        label = labels_unique
      )

    ### NULLIFY for CRAN check
    center <- label <- cut <- NULL

    if (cluster_labels) {
      # Place cluster labels along diagonal
      # if (!flip) {
      #   g <- g + ggplot2::geom_label(data = clusters,
      #     ggplot2::aes(
      #       x = center,
      #       y = nrow(m) - center,
      #       label = label
      #     ))
      # } else{
      #   g <- g + ggplot2::geom_label(data = clusters,
      #     ggplot2::aes(
      #       x = ncol(m) - center,
      #       y = nrow(m) - center,
      #       label = label
      #     ))
      # }

      # Place cluster labels on top as x-axis (needs the colnames set as a sequence)
      if (reverse_columns) {
        suppressMessages(
          g <-
            g + ggplot2::scale_x_discrete(
              breaks = ncol(m) - clusters$center,
              label = clusters$label,
              expand = c(0, 0),
              position = "top"
            ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(
                angle = 0,
                vjust = 0.5,
                hjust = .5
              )
            ) +
            ggplot2::labs(x = "Cluster")
        )
      } else{
        suppressMessages(
          g <- g + ggplot2::scale_x_discrete(
            breaks = clusters$center,
            label = clusters$label,
            expand = c(0, 0),
            position = "top"
          ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(
                angle = 0,
                vjust = 0.5,
                hjust = .5
              )
            ) +
            ggplot2::labs(x = "Cluster")
        )
      }

      if (cluster_lines) {
        ## draw lines separating the clusters

        if (reverse_columns) {
          g <-
            g + ggplot2::geom_hline(data = clusters, ggplot2::aes(yintercept = nrow(m) - cut + .5)) +
            ggplot2::geom_vline(data = clusters, ggplot2::aes(xintercept = ncol(m) - cut + .5))
        } else{
          g <-
            g + ggplot2::geom_hline(data = clusters, ggplot2::aes(yintercept = nrow(m) - cut + .5)) +
            ggplot2::geom_vline(data = clusters, ggplot2::aes(xintercept = cut + .5))
        }
      }
    }
  }

  # reverse color
  suppressMessages(g <-
      g + .gg_sequential_pal(dist = TRUE))

  g
}
