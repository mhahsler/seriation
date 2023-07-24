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


## Hierarchical clustering related seriations
.hc_control <- list(hclust = NULL,
                    linkage = "complete")

attr(.hc_control, "help") <-
  list(hclust = "a precomputed hclust object (optional)",
       linkage = "hclust method")

.hclust_helper <- function(d, control = NULL) {
  # Deprecated method control argument
  if (!is.null(control$method)) {
    warning("control parameter method is deprecated. Use linkage instead!")
    control$linkage <- control$method
    control$method <- NULL
  }

  control <- .get_parameters(control, .hc_control)

  if (!is.null(control$hclust))
    return(control$hclust)
  return(hclust(d, method = control$linkage))
}

seriate_dist_hc <- function(x, control = NULL)
  .hclust_helper(x, control)
seriate_dist_hc_single <- function(x, control = NULL)
  .hclust_helper(x, control = list(linkage = "single"))
seriate_dist_hc_average <- function(x, control = NULL)
  .hclust_helper(x, control = list(linkage = "average"))
seriate_dist_hc_complete <- function(x, control = NULL)
  .hclust_helper(x, control = list(linkage = "complete"))
seriate_dist_hc_ward <- function(x, control = NULL)
  .hclust_helper(x, control = list(linkage = "ward.D2"))

seriate_dist_gw <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method = "GW")
seriate_dist_gw_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method = "GW")
seriate_dist_gw_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method = "GW")
seriate_dist_gw_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method = "GW")
seriate_dist_gw_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method = "GW")

seriate_dist_olo <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method = "OLO")
seriate_dist_olo_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method = "OLO")
seriate_dist_olo_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method = "OLO")
seriate_dist_olo_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method = "OLO")
seriate_dist_olo_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method = "OLO")


.hc_desc <-
  "Using the order of the leaf nodes in a dendrogram obtained by hierarchical clustering"

.optHCPL <- .opt("Path_length", "restricted by dendrogram")

set_seriation_method("dist", "HC", seriate_dist_hc,
                     .hc_desc, .hc_control)

set_seriation_method("dist",
                     "HC_single",
                     seriate_dist_hc_single,
                     paste(.hc_desc, "(single link)"))

set_seriation_method(
  "dist",
  "HC_complete",
  seriate_dist_hc_complete,
  paste(.hc_desc, "(complete link).")
)

set_seriation_method("dist",
                     "HC_average",
                     seriate_dist_hc_average,
                     paste(.hc_desc, "(avg. link)."))

set_seriation_method("dist",
                     "HC_ward",
                     seriate_dist_hc_ward,
                     paste(.hc_desc, "(Ward's method)."))

.gw_desc <-
  "Using the order of the leaf nodes in a dendrogram obtained by hierarchical clustering and reordered by the Gruvaeus and Wainer (1972) heuristic"

set_seriation_method("dist",
                     "GW",
                     seriate_dist_gw,
                     .gw_desc,
                     .hc_control,
                     optimizes = .optHCPL)

set_seriation_method(
  "dist",
  "GW_single",
  seriate_dist_gw_single,
  paste(.gw_desc, "(single link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "GW_average",
  seriate_dist_gw_average,
  paste(.gw_desc, "(avg.link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "GW_complete",
  seriate_dist_gw_complete,
  paste(.gw_desc, "(complete link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "GW_ward",
  seriate_dist_gw_ward,
  paste(.gw_desc, "(Ward's method)"),
  optimizes = .optHCPL
)

.olo_desc <-
  "Using the order of the leaf nodes in a dendrogram obtained by hierarchical clustering and reordered by with optimal leaf ordering (Bar-Joseph et al., 2001)"

set_seriation_method("dist",
                     "OLO",
                     seriate_dist_olo,
                     .olo_desc,
                     .hc_control,
                     optimizes = .optHCPL)

set_seriation_method(
  "dist",
  "OLO_single",
  seriate_dist_olo_single,
  paste(.olo_desc, "(single link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "OLO_average",
  seriate_dist_olo_average,
  paste(.olo_desc, "(avg. link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "OLO_complete",
  seriate_dist_olo_complete,
  paste(.olo_desc, "(complete link)"),
  optimizes = .optHCPL
)

set_seriation_method(
  "dist",
  "OLO_ward",
  seriate_dist_olo_ward,
  paste(.olo_desc, "(Ward's method)"),
  optimizes = .optHCPL
)
