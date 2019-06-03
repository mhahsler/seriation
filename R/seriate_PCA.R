#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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


## use the projection on the first pricipal component to determine the
## order
seriate_matrix_fpc <- function(x, control = NULL) {
  control <- .get_parameters(control, list(
    center = TRUE,
    scale. = FALSE,
    tol = NULL,
    verbose = FALSE
  ))

  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol
  verbose <- control$verbose

  pr <- prcomp(x, center = center, scale. = scale., tol = tol)
  scores <- pr$x[,1]
  row <- order(scores)
  if(verbose) cat("row: first principal component explains",
    pr$sdev[1] / sum(pr$sdev)* 100,"%\n")

  pr <- prcomp(t(x), center = center, scale. = scale., tol = tol)
  scores <- pr$x[,1]
  col <- order(scores)
  if(verbose) cat("col: first principal component explains",
    pr$sdev[1] / sum(pr$sdev)* 100,"%\n")

  names(row) <- rownames(x)[row]
  names(col) <- colnames(x)[col]

  list(row = row, col = col)
}

## Angle between the first 2 PCS. Fiendly (2002)
.order_angle <- function(x) {

  alpha <- atan2(x[,1], x[,2])
  o <- order(alpha)
  cut <- which.max(abs(diff(c(alpha[o], alpha[o[1]]+2*pi))))
  if(cut==length(o)) o
  else o[c((cut+1):length(o), 1:(cut))]

}


seriate_matrix_angle <- function(x, control = NULL) {
  control <- .get_parameters(control, list(
    center = TRUE,
    scale. = FALSE,
    tol = NULL
  ))

  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol

  pr <- prcomp(x, center = center, scale. = scale., tol = tol)
  row <- .order_angle(pr$x[,1:2])

  pr <- prcomp(t(x), center = center, scale. = scale., tol = tol)
  col <- .order_angle(pr$x[,1:2])

  names(row) <- rownames(x)[row]
  names(col) <- colnames(x)[col]

  list(row = row, col = col)
}

set_seriation_method("matrix", "PCA", seriate_matrix_fpc,
  "First principal component")
set_seriation_method("matrix", "PCA_angle", seriate_matrix_angle,
  "First two principal components (angle)")

