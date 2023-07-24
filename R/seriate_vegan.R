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

.monoMDS_control <- structure({
  l <- as.list(args(vegan::monoMDS))
  l$k <- NULL
  l$model <- "global"
  tail(head(l,-2L),-1L)
}, help = list(y = "See ? monoMDS for help"))

seriate_dist_monoMDS <- function(x, control = NULL) {
  control <- .get_parameters(control, .monoMDS_control)
  r <- do.call(vegan::monoMDS, c(list(x, k = 1), control))
  conf <- r$points

  if (control$verbose) {
    r$call <- NULL
    print(r)
  }

  structure(order(conf), configutation = conf)
}

set_seriation_method(
  "dist",
  "monoMDS",
  seriate_dist_monoMDS,
  "Kruskal's (1964a,b) non-metric multidimensional scaling (NMDS) using monotone regression.",
  control = .monoMDS_control,
  randomized = TRUE,
  optimizes = .opt("MDS_stress", "Kruskal's monotone regression stress")
)



.isomap_control <- structure(
  list(k = 30,
       path = "shortest"),
  help = list(k = "number of shortest dissimilarities retained for a point",
              path = "method used in to estimate the shortest path (\"shortest\"/\"extended\")")
)

seriate_dist_isomap <- function(x, control = NULL) {
  control <- .get_parameters(control, .isomap_control)
  r <- do.call(vegan::isomap, c(list(x, ndim = 1), control))
  conf <- r$points

  if (control$verbose) {
    r$call <- NULL
    print(r)
  }

  structure(order(conf), configutation = conf)
}

set_seriation_method(
  "dist",
  "isomap",
  seriate_dist_isomap,
  "Isometric feature mapping ordination",
  control = .isomap_control,
  optimizes = .opt(NA, "Stress on shortest path distances")
)

.metaMDS_control <- structure({
  l <- as.list(args(vegan::metaMDS))
  l <- tail(head(l, -2L), -1L)
  l$k <- NULL
  l$engine <- "monoMDS"
  l$noshare <- FALSE
  #l$distance = "euclidean"
  l$trace <- 0
  l$verbose <- FALSE
  l
  }, help = list(distance = "see ? metaMDS for help")
)

seriate_dist_metaMDS <- function(x, control = NULL) {
  control <- .get_parameters(control, .metaMDS_control)

  r <- do.call(vegan::metaMDS, c(list(x, k = 1), control))
  conf <- r$points

  if(control$verbose && control$trace == 0)
    control$trace <- 1

  if (control$verbose) {
    r$call <- NULL
    r$data <- NULL
    print(r)
  }

  structure(order(conf), configutation = conf)
}

set_seriation_method(
  "dist",
  "metaMDS",
  seriate_dist_metaMDS,
  "Nonmetric Multidimensional Scaling with Stable Solution from Random Starts.",
  control = .metaMDS_control,
  randomized = FALSE,          ### it is randomized, but internally does replication
  optimizes = .opt("MDS_stress", "Kruskal's monotone regression stress")
)
