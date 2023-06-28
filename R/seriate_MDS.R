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


.mds_control <- list(type = "cmdscale",
                     tol = 1e-3,
                     maxit = 100
)

seriate_dist_mds <- function(x, control = NULL) {
  ### accept deprecated method
  if (!is.null(control$method)) {
    control$type <- control$method
    control$method <- control$type
    warning("seriation method mds: control parameter method is deprecated. Use type instead!")
  }

  control <- .get_parameters(control, .mds_control)

  method <-
    match.arg(control$type, c("cmdscale", "isoMDS", "sammon"))

  if (method == "cmdscale") {
    sc <- stats::cmdscale(x, k = 1)
    o <- order(sc[, 1])
    attr(o, "embedding") <- sc[, 1]
    return(o)

  } else if (method == "isoMDS") {
    sc <-
      MASS::isoMDS(
        x + 1e-6,
        trace = FALSE,
        k = 1,
        maxit = control$maxit,
        tol = control$tol
      )
    o <- order(sc$points[, 1])
    attr(o, "embedding") <- sc$points[, 1]
    return(o)

  } else if (method == "sammon") {
    sc <- MASS::sammon(
      x + 1e-6,
      trace = FALSE,
      k = 1,
      niter = control$maxit,
      tol = control$tol
    )
    o <- order(sc$points[, 1])
    attr(o, "embedding") <- sc$points[, 1]
    return(o)

  } else
    stop("unknown method.")

}

seriate_dist_mds_metric <- function(x, control = NULL)
  seriate_dist_mds(x, control = list(type = "cmdscale"))

seriate_dist_mds_nonmetric <- function(x, control = NULL)
  seriate_dist_mds(x, control = list(type = "isoMDS"))

seriate_dist_mds_sammon <- function(x, control = NULL)
  seriate_dist_mds(x, control = list(type = "sammon"))

## Angle between the first 2 PCS. Friendly (2002)
seriate_dist_angle <- function(x, control = NULL) {
  .get_parameters(control, NULL)

  sc <- stats::cmdscale(x, k = 2)
  .order_angle(sc)
}


set_seriation_method(
  "dist",
  "MDS",
  seriate_dist_mds,
  "Order using the first component found by multidimensional scaling. Element type in control can be \"cmdscale\", \"isoMDS\" or \"sammon\".",
  .mds_control,
  optimizes = "Stress (Moore stress)"
)

set_seriation_method(
  "dist",
  "MDS_metric",
  seriate_dist_mds_metric,
  "Order using the first component found by classical metric multidimensional scaling (cmdscsale).",
  optimizes = "Stress (Moore stress)"
)

set_seriation_method(
  "dist",
  "MDS_nonmetric",
  seriate_dist_mds_nonmetric,
  "Order using the first component found by Kruskal's non-metric multidimensional scaling (isoMDS).",
  optimizes = "Stress (Moore stress)"
)

set_seriation_method(
  "dist",
  "MDS_sammon",
  seriate_dist_mds_sammon,
  "Order using the first component found by Sammon's non-linear mapping (sammon).",
  optimizes = "Stress (Moore stress)"
)

set_seriation_method(
  "dist",
  "MDS_angle",
  seriate_dist_angle,
  "Order by the angle in this space given by the first two components found by metric MDS (Friendly, 2002)."
)
