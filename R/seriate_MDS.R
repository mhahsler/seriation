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


.mds_control <- list(
  method = "cmdscale"
)

## Multidimensional scaling
seriate_dist_mds <- function(x, control = NULL){
  control <- .get_parameters(control, .mds_control)

  if(control$method == "cmdscale" ) {
    sc <- cmdscale(x, k=1)
    return(order(sc[,1]))

  }else if(control$method == "isoMDS"){
    sc <- MASS::isoMDS(x+1e-6, trace = FALSE, k=1)
    return(order(sc$points[,1]))

  }else if(control$method == "sammon") {
    sc <- MASS::sammon(x+1e-6, trace = FALSE, k=1)
    return(order(sc$points[,1]))

  }else stop("unknown method")

}

seriate_dist_mds_metric <- function(x, control = NULL)
  seriate_dist_mds(x, control=list(method="cmdscale"))

seriate_dist_mds_nonmetric <- function(x, control = NULL)
  seriate_dist_mds(x, control=list(method="isoMDS"))

## Angle between the first 2 PCS. Fiendly (2002)
seriate_dist_angle <- function(x, control = NULL) {
  .get_parameters(control, NULL)

  sc <- cmdscale(x, k=2)
  .order_angle(sc)
}


set_seriation_method("dist", "MDS", seriate_dist_mds,
  "Order using the first component found by multidimensional scaling. Element method in control can be \"cmdscale\", \"isoMDS\" or \"sammon\".", .mds_control)
set_seriation_method("dist", "MDS_metric", seriate_dist_mds_metric,
  "Order using the first component found by metric multidimensional scaling (cmdscsale).")
set_seriation_method("dist", "MDS_nonmetric", seriate_dist_mds_nonmetric,
  "Order using the first component found by non-metric multidimensional scaling (isoMDS).")
set_seriation_method("dist", "MDS_angle", seriate_dist_angle,
  "Order by the angle in this space given by the first two components found by MDS (Friendly, 2002).")

