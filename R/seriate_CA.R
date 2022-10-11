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

## use the projection on the first principal component to determine the
## order
.ca_contr <- list(
  dim = 1L,
  ca_param = NULL
)

seriate_matrix_ca <- function(x, control = NULL) {
  control <- .get_parameters(control, .ca_contr)

  mat.ca <- do.call(ca::ca, c(list(obj = x), control$ca_param))
  rcoord <- mat.ca$rowcoord    # row coordinates
  row <- order(rcoord[, control$dim])
  ccoord <- mat.ca$colcoord    # col coordinates
  col <- order(ccoord[, control$dim])

  #names(row) <- rownames(x)[row]
  #names(col) <- colnames(x)[col]

  list(row = row, col = col)
}


set_seriation_method(
  "matrix",
  "CA",
  seriate_matrix_ca,
  "This method calculates a correspondence analysis of the matrix and computes an order according to the scores on a correspondence analysis dimension.",
  .ca_contr
)
