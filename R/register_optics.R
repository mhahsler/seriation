#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
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

## register seriation based on OPTICS

register_optics <- function() {
  check_installed("dbscan")

  .contr <- list(
    eps = NULL,
    minPts = 5
  )

  optics_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    dbscan::optics(x, eps = control$eps, minPts = control$minPts)$order
  }

  set_seriation_method(
    "dist",
    "optics",
    optics_order,
    "Use ordering points to identify the clustering structure (OPTICS) to create an order",
    .contr
  )
}
