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


seriate_dist_random <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  o <- 1:attr(x, "Size")
  sample(o)
}

seriate_matrix_random <-
  function(x, control, margin = seq_along(dim(x))) {
    control <- .get_parameters(control, NULL)
    lapply(seq_along(dim(x)), function(i)
      if (i %in% margin)
        sample(seq(dim(x)[i]))
      else
        NA)
  }

set_seriation_method("dist",
                     "Random",
                     seriate_dist_random,
                     "Random permutation",
                     randomized = TRUE,
                     optimized = "None")

set_seriation_method("matrix",
                     "Random",
                     seriate_matrix_random,
                     "Random permutation",
                     randomized = TRUE,
                     optimized = "None")

set_seriation_method("array",
                     "Random",
                     seriate_matrix_random,
                     "Random permutation",
                     randomized = TRUE,
                     optimized = "None")
