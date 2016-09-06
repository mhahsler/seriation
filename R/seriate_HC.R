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


## Hierarchical clustering related seriations
.hclust_helper <- function(d, control = NULL){
  control <- .get_parameters(control, list(
    hclust = NULL,
    method = "average"
    ))

  if(!is.null(control$hclust)) return(control$hclust)
  return(hclust(d, method = control$method))
}

seriate_dist_hc <- function(x, control = NULL) .hclust_helper(x, control)
seriate_dist_hc_single <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="single"))
seriate_dist_hc_average <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="average"))
seriate_dist_hc_complete <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="complete"))
seriate_dist_hc_ward <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="ward.D2"))

## workhorses are in seriation.hclust
seriate_dist_gw <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method="GW")
seriate_dist_gw_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method="GW")
seriate_dist_gw_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method="GW")
seriate_dist_gw_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method="GW")
seriate_dist_gw_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method="GW")


seriate_dist_olo <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method="OLO")
seriate_dist_olo_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method="OLO")
seriate_dist_olo_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method="OLO")
seriate_dist_olo_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method="OLO")
seriate_dist_olo_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method="OLO")

set_seriation_method("dist", "HC", seriate_dist_hc,
  "Hierarchical clustering")
set_seriation_method("dist", "HC_single", seriate_dist_hc_single,
  "Hierarchical clustering (single link)")
set_seriation_method("dist", "HC_complete", seriate_dist_hc_complete,
  "Hierarchical clustering (complete link)")
set_seriation_method("dist", "HC_average", seriate_dist_hc_average,
  "Hierarchical clustering (avg. link)")
set_seriation_method("dist", "HC_ward", seriate_dist_hc_ward,
  "Hierarchical clustering (Ward's method)")

set_seriation_method("dist", "GW", seriate_dist_gw,
  "Hierarchical clustering reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_single", seriate_dist_gw_single,
  "Hierarchical clustering (single link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_average", seriate_dist_gw_average,
  "Hierarchical clustering (avg. link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_complete", seriate_dist_gw_complete,
  "Hierarchical clustering (complete link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_ward", seriate_dist_gw_ward,
  "Hierarchical clustering (Ward's method) reordered by Gruvaeus and Wainer heuristic")


set_seriation_method("dist", "OLO", seriate_dist_olo,
  "Hierarchical clustering (single link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_single", seriate_dist_olo_single,
  "Hierarchical clustering with optimal leaf ordering")
set_seriation_method("dist", "OLO_average", seriate_dist_olo_average,
  "Hierarchical clustering (avg. link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_complete", seriate_dist_olo_complete,
  "Hierarchical clustering (complete link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_ward", seriate_dist_olo_ward,
  "Hierarchical clustering (Ward's method) with optimal leaf ordering")

