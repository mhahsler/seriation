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


## use the projection on the first principal component to determine the
## order
.pca_contr <- list(
  center = TRUE,
  scale = FALSE,
  verbose = FALSE
)
attr(.pca_contr, "help") <- list(
  center = "center the data (mean = 0)?",
  scale = "scale to unit variance?",
  verbose = FALSE
)

seriate_matrix_fpc <- function(x, control = NULL, margin) {
  control <- .get_parameters(control, .pca_contr)

  center  <- control$center
  scale  <- control$scale
  verbose <- control$verbose

  o <- list(row = NA, col = NA)

  if (1L %in% margin) {
    pr <- stats::prcomp(x,
                        center = center,
                        scale. = scale,
                        rank. = 1L)
    scores <- pr$x[, 1]
    os <- order(scores)
    o$row <- structure(os, names = rownames(x)[os], configuration = scores)

    if (verbose)
      cat("Rows: first PC explains",
          pr$sdev[1] / sum(pr$sdev) * 100,
          "%\n")
  }

  if (2L %in% margin) {
    x <- t(x)
    pr <- stats::prcomp(x,
                 center = center,
                 scale. = scale,
                 rank. = 1L)
    scores <- pr$x[, 1]
    os <- order(scores)
    o$col <- structure(os, names = rownames(x)[os], configuration = scores)

     if (verbose)
      cat("Cols: first PC explains",
          pr$sdev[1] / sum(pr$sdev) * 100,
          "%\n")
  }

  if (verbose)
    cat("\n")

  o
}


seriate_matrix_angle <- function(x, control = NULL, margin) {
  control <- .get_parameters(control, .pca_contr)

  center  <- control$center
  scale  <- control$scale

  if (nrow(x) < 2L || ncol(x) < 2L)
    stop("PCA angle needs at least 2 rows and 2 columns!")

  if (1L %in% margin) {
      pr <- prcomp(x,
                   center = center,
                   scale. = scale,
                   rank = 2L)
      row <- .order_angle(pr$x[, 1:2])
  } else
    row <- NA

  if (2L %in% margin) {
      pr <- prcomp(t(x),
                   center = center,
                   scale. = scale,
                   rank = 2L)
      col <- .order_angle(pr$x[, 1:2])
  } else
    col <- NA

  list(row = row, col = col)
}

set_seriation_method(
  "matrix",
  "PCA",
  seriate_matrix_fpc,
  "Uses the projection of the data on its first principal component to determine the order.",
  .pca_contr,
  optimizes = .opt(NA, "Least squares for each dimension (for Euclidean distances).")
)

set_seriation_method(
  "matrix",
  "PCA_angle",
  seriate_matrix_angle,
  "Uses the angular order in the 2D PCA projection space split by the larges gap.",
  .pca_contr
)
