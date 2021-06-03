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



## Cluster visualization by proximity matrix shading

ggdissplot <- function(x,
  labels = NULL,
  method = "Spectral",
  control = NULL,
  options = NULL,
  ...) {
  check_installed("ggplot2")

  # add ... to options
  options <- c(options, list(...))

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

  options <- .get_parameters(
    options,
    list(
      cluster_labels = TRUE,
      lines       = TRUE,
      averages	  = c(FALSE, TRUE),
      labRow = NULL,
      labCol = NULL,
      flip		    = FALSE,
      # silhouettes = FALSE, ???
      threshold   = NULL
    )
  )

  m       <- as.matrix(x$x_reordered)
  k       <- x$k
  dim     <- attr(x$x_reordered, "Size")
  labels  <- x$labels
  labels_unique <- unique(labels)

  if (is.na(options$averages[1]))
    m[upper.tri(m)] <- NA
  if (is.na(options$averages[2]))
    m[lower.tri(m)] <- NA
  options$averages[is.na(options$averages)] <- FALSE

  if (!is.null(x$cluster_dissimilarities)
    && !is.null(labels)
    && any(options$averages)) {
    for (i in 1:k) {
      for (j in 1:k) {
        ## check empty clusters
        if (is.na(labels_unique[i]))
          next
        if (is.na(labels_unique[j]))
          next

        ## upper panels stay unchanged
        if (i < j && options$averages[1]) {
          m[labels == labels_unique[i], labels == labels_unique[j]] <-
            x$cluster_dissimilarities[i, j]
        }

        ## do lower panels
        if (i > j && options$averages[2]) {
          m[labels == labels_unique[i], labels == labels_unique[j]] <-
            x$cluster_dissimilarities[i, j]
        }

        ## do diagonal
        if (i == j) {
          block <- m[labels == labels_unique[i],
            labels == labels_unique[j]]


          if (options$averages[1]) {
            block[upper.tri(block, diag = TRUE)] <-
              x$cluster_dissimilarities[i, j]

            m[labels == labels_unique[i],
              labels == labels_unique[j]] <- block
          }

          if (options$averages[2]) {
            block[lower.tri(block, diag = TRUE)] <-
              x$cluster_dissimilarities[i, j]

            m[labels == labels_unique[i],
              labels == labels_unique[j]] <- block
          }

        }
      }
    }
  }

  if (options$flip)
    m <- m[, ncol(m):1]

  g <- ggpimage(m, labRow = options$labRow, labCol = options$labCol)

  # add cluster lines and labels
  if (!is.null(labels)) {
    cluster_width		<- tabulate(labels)[labels_unique]
    cluster_cuts		<- cumsum(cluster_width)
    cluster_center	<- cluster_cuts - cluster_width / 2

    clusters <-
      data.frame(
        center = cluster_center,
        cut = cluster_cuts,
        with = cluster_width,
        label = labels_unique
      )

    ### NULLIFY for CRAN check
    center <- label <- cut <- NULL

    if (options$cluster_labels) {
      ## TODO: above the plot
      if (!options$flip) {
        g <- g + ggplot2::geom_text(data = clusters,
          ggplot2::aes(
            x = center,
            y = nrow(m) - center,
            label = label
          ))
      } else{
        g <- g + ggplot2::geom_text(data = clusters,
          ggplot2::aes(
            x = ncol(m) - center,
            y = nrow(m) - center,
            label = label
          ))
      }

      if (options$lines) {
        ## draw lines separating the clusters

        if (!options$flip) {
          g <-
            g + ggplot2::geom_hline(data = clusters, ggplot2::aes(yintercept = nrow(m) - cut + .5)) +
            ggplot2::geom_vline(data = clusters, ggplot2::aes(xintercept = cut + .5))
        } else{
          g <-
            g + ggplot2::geom_hline(data = clusters, ggplot2::aes(yintercept = nrow(m) - cut + .5)) +
            ggplot2::geom_vline(data = clusters, ggplot2::aes(xintercept = ncol(m) - cut + .5))
        }
      }
    }
  }
  g

}
