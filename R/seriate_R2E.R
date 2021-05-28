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

## uses a sequence of correlation matrices and finds  the first matrix
## with rank 2. The elements are projected into the plane spanned by the
## first two eigenvectors. All points are lying on a ellipse. The order
## of the elements on the ellipse is returned (see Chen 2002).
seriate_dist_chen <- function(x, control = NULL) {
  .get_parameters(control, NULL)

  x <- as.matrix(x)

  rank <- qr(x)$rank

  ## find the first correlation matrix of rank 2
  n <- 0
  while (rank > 2) {
    x <- cor(x)
    n <- n + 1
    rank <- qr(x)$rank
  }

  ## project the matrix on the first 2 eigenvectors
  e <- eigen(x)$vectors[, 1:2]

  ## extract the order
  ## chen says that he uses the one of the two possible cuts
  ## that separate the points at rank 1. Since the points just
  ## separate further towards right and left, cutting on the vertical
  ## axis of the ellipse yields the same result.

  right <- which(e[, 1] >= 0)
  right <- right[order(e[right, 2], decreasing = TRUE)]
  left <- which(e[, 1] < 0)
  left <- left[order(e[left, 2])]

  o <- c(right, left)
  names(o) <- labels(x)[o]
  o
}

#set_seriation_method("dist", "Chen", seriate_dist_chen,
#  "Rank-two ellipse seriation")

set_seriation_method("dist", "R2E", seriate_dist_chen,
  "Rank-two ellipse seriation")
