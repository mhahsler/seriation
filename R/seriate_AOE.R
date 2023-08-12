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

is_correlation_matrix <- function(x) {
  if(!isSymmetric(x))
    return (FALSE)
  if(any(diag(x) != 1))
    return (FALSE)
  if(any(x > 1))
    return (FALSE)
  if(any(x < -1))
    return (FALSE)

  return(TRUE)
}

# AOE for correlation matrices
seriate_corr_matrix_AOE <- function(x, control = NULL, margin) {
  if(!is_correlation_matrix(x)) {
    warning("x is not a correlation matrix. Using method 'PCA_angle' instead.")
    return(seriate_matrix_angle(x, control, margin))
  }

  sc <- eigen(x)$vectors[, 1:2]
  o <- .order_angle(sc)

  list(row = o, col = o)
}

## Angle between the first 2 PCs.
# Friendly, M. (2002), "Corrgrams: Exploratory Displays for Correlation Matrices," The American Statistician,56, 316-324.
# Friendly, M. and Kwan, E. (2003), "Effect ordering for data displays," Computational Statistics & Data Analysis, 43, 509-539.
.order_angle <- function(x) {
  alpha <- atan2(x[, 1], x[, 2])
  o <- order(alpha)

  # cut at largest gap. Gap is before!
  cut <- which.max(abs(diff(c(
    alpha[o], alpha[o[1]] + 2 * pi, lag = 1
  ))))

  if (cut != 1L && cut != length(o) + 1L)
    o <- o[c((cut):length(o), 1:(cut - 1L))]
  o

}

set_seriation_method(
  "matrix",
  "AOE",
  seriate_corr_matrix_AOE,
  "Order by the angle of the first two eigenvectors (for correlation matrices)",
)

