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
  scale. = FALSE,
  tol = NULL,
  verbose = FALSE
)

seriate_matrix_fpc <- function(x, control = NULL, margin = NULL) {
  control <- .get_parameters(control, .pca_contr)

  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol
  verbose <- control$verbose

  if (1L %in% margin) {
    pr <- stats::prcomp(x,
                        center = center,
                        scale. = scale.,
                        tol = tol)
    scores <- pr$x[, 1]
    row <- order(scores)
    attr(row, "embedding") <- scores

    if (verbose)
      cat("Rows: first PC explains",
          pr$sdev[1] / sum(pr$sdev) * 100,
          "%\n")
    names(row) <- rownames(x)[row]
  } else
    row <- NA


  if (2L %in% margin) {
    pr <- stats::prcomp(t(x),
                 center = center,
                 scale. = scale.,
                 tol = tol)
    scores <- pr$x[, 1]
    col <- order(scores)
    attr(col, "embedding") <- scores

     if (verbose)
      cat("Cols: first PC explains",
          pr$sdev[1] / sum(pr$sdev) * 100,
          "%\n")
    names(col) <- colnames(x)[col]
  } else
    col <- NA

  if (verbose)
    cat("\n")

  list(row = row, col = col)
}

## Angle between the first 2 PCs.
# Friendly, M. (2002), "Corrgrams: Exploratory Displays for Correlation Matrices," The American Statistician,56, 316-324.
# Friendly, M. and Kwan, E. (2003), "Eect ordering for data displays," Computational Statistics & Data Analysis, 43, 509-539.
.order_angle <- function(x) {
  alpha <- atan2(x[, 1], x[, 2])
  o <- order(alpha)
  cut <- which.max(abs(diff(c(
    alpha[o], alpha[o[1]] + 2 * pi
  ))))
  if (cut == length(o))
    o
  else
    o[c((cut + 1):length(o), 1:(cut))]

}

.angle_contr <- list(center = TRUE,
                     scale. = FALSE,
                     tol = NULL)

seriate_matrix_angle <- function(x, control = NULL, margin = seq_along(dim(x))) {
  control <- .get_parameters(control, .angle_contr)

  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol

  if (1L %in% margin) {
    pr <- prcomp(x,
                 center = center,
                 scale. = scale.,
                 tol = tol)
    row <- .order_angle(pr$x[, 1:2])
    names(row) <- rownames(x)[row]
  } else
    row <- NA

  if (2L %in% margin) {
    pr <- prcomp(t(x),
                 center = center,
                 scale. = scale.,
                 tol = tol)
    col <- .order_angle(pr$x[, 1:2])
    names(col) <- colnames(x)[col]
  } else
    col <- NA


  list(row = row, col = col)
}

set_seriation_method(
  "matrix",
  "PCA",
  seriate_matrix_fpc,
  "Uses the projection of the data on its first principal component to determine the order. Note that for a distance matrix calculated from x with Euclidean distance, this methods minimizes the least square criterion.",
  .pca_contr,
  optimizes = "Least squares for each dimension."
)

set_seriation_method(
  "matrix",
  "PCA_angle",
  seriate_matrix_angle,
  "Projects the data on the first two principal components and then orders by the angle in this space. The order is split by the larges gap between adjacent angles.",
  .angle_contr
)
