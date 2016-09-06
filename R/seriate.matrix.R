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



## seriate matrices

seriate.matrix <- function(x, method = "PCA", control = NULL,
  margin = c(1,2), ...)
  .seriate_array_helper(x, method, control, margin,
    datatype = "matrix", defmethod = "BEA_TSP", ...)


seriate_matrix_identity <- function(x, control) {
  control <- .get_parameters(control, NULL)

  l <- lapply(dim(x), seq)
  for(i in 1:length(dim(x))) names(l[[i]]) <- labels(x)[[i]]
  l
}

seriate_matrix_random <- function(x, control) {
  control <- .get_parameters(control, NULL)

  l <- lapply(dim(x), FUN = function(l) sample(seq(l)))
  for(i in 1:length(dim(x))) names(l[[i]]) <- labels(x)[[i]][l[[i]]]
  l
}


set_seriation_method("matrix", "Identity", seriate_matrix_identity,
  "Identity permutation")
set_seriation_method("matrix", "Random", seriate_matrix_random,
  "Random permutation")


## these also work for general arrays!
set_seriation_method("array", "Identity", seriate_matrix_identity,
  "Identity permutation")
set_seriation_method("array", "Random", seriate_matrix_random,
  "Random permutation")

