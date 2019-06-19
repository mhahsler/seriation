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

## Bridge to package tsp
.tsp_control <- list(
  method="arbitrary insertion",
  rep = 10,
  two_opt = TRUE
)

seriate_dist_tsp <- function(x, control = NULL){
  ## add a dummy city for cutting
  tsp <- insert_dummy(TSP(x), n = 1, label = "cut_here")

  if(is.null(control)) control <- .tsp_control

  tour <- solve_TSP(tsp, method = control$method,
    control = control)

  o <- cut_tour(tour, cut = "cut_here", exclude_cut = TRUE)
  names(o) <- labels(x)[o]
  o
}

set_seriation_method("dist", "TSP", seriate_dist_tsp,
  "Minimize Hamiltonian path length with a TSP solver (see solve_TSP in package TSP for available methods).", .tsp_control)

