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

# controls are only used for isoMDS and sammon
.mds_control <- list(add = FALSE)

seriate_dist_mds <- function(x, control = NULL) {
  ### accept deprecated method
  if (!is.null(control$method)) {
    control$method <- NULL
    warning("seriation method mds: control parameter method is deprecated and ignored!")
  }

  control <- .get_parameters(control, .mds_control)

  # eig = TRUE makes sure we get a list back
  sc <- stats::cmdscale(x, k = 1,  eig = TRUE, add = control$add)
  sc <- drop(sc$points)
  o <- order(sc)
  attr(o, "configuration") <- sc
  o
}


.mds_isoMDS_control <- list(
  add = 1e-9,
  # to avoid 0 distances
  maxit = 50,
  trace = FALSE,
  tol = 1e-3,
  p = 2
)

seriate_dist_mds_isoMDS <- function(x, control = NULL) {
  control <- .get_parameters(control, .mds_isoMDS_control)

  sc <-
    MASS::isoMDS(
      x + control$add,
      k = 1,
      maxit = control$maxit,
      trace = control$trace,
      tol = control$tol,
      p = control$p
    )
  o <- order(sc$points[, 1])
  attr(o, "configuration") <- sc$points[, 1]
  o
}

.mds_sammon_control <- list(
  add = 1e-9,
  # to avoid 0 distances
  niter = 100,
  trace = FALSE,
  magic = 0.2,
  tol = 1e-4
)

seriate_dist_mds_sammon <- function(x, control = NULL) {
  control <- .get_parameters(control, .mds_sammon_control)

  sc <- MASS::sammon(
    x + control$add,
    y = jitter(stats::cmdscale(x, k = 1)),
    ### fixes issue with duplicates
    k = 1,
    niter = control$niter,
    trace = control$trace,
    magic = control$magic,
    tol = control$tol
  )

  o <- order(sc$points[, 1])
  attr(o, "configuration") <- sc$points[, 1]
  o
}

## Angle between the first 2 PCS. Friendly (2002)
seriate_dist_angle <- function(x, control = NULL) {
  control <- .get_parameters(control, .mds_control)

  sc <- stats::cmdscale(x, k = 2, eig = TRUE, add = control$add)
  sc <- sc$points
  o <- .order_angle(sc)
  attr(o, "configuration") <- sc
  o
}


set_seriation_method(
  "dist",
  "MDS",
  seriate_dist_mds,
  "Order along the 1D classical metric multidimensional scaling",
  control = .mds_control,
  optimizes = "Other (MDS strain)"
)

set_seriation_method(
  "dist",
  "MDS_isoMDS",
  seriate_dist_mds_isoMDS,
  "Order along the 1D Kruskal's non-metric multidimensional scaling",
  control = .mds_isoMDS_control,
  optimizes = "Other (MDS stress with monotonic transformation)"
)

set_seriation_method(
  "dist",
  "MDS_sammon",
  seriate_dist_mds_sammon,
  "Order along the 1D Sammon's non-linear mapping",
  control = .mds_sammon_control,
  optimizes = "Other (scale free, weighted MDS stress)"
)

set_seriation_method(
  "dist",
  "MDS_angle",
  seriate_dist_angle,
  "Order by the angle in this space given by 2D metric MDS.",
  control = .mds_control
)
