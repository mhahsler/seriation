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


#' Register Seriation Based on OPTICS
#'
#' Use ordering points to identify the clustering structure (OPTICS) for [seriate()].
#'
#' Registers the method \code{"optics"} for [seriate()]. This method applies
#' the OPTICS ordering algorithm to create an ordering.
#'
#' \bold{Note:} Package \pkg{dbscan} needs to be installed.
#'
#' @aliases register_optics optics OPTICS
#' @seealso [dbscan::optics()] in \pkg{dbscan}.
#' @family seriation
#' @returns Nothing.
#'
#' @references Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg
#' Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure.
#' ACM SIGMOD international conference on Management of data. ACM Press. pp.
#' 49-60. \doi{10.1145/304181.304187}
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_optics()
#' get_seriation_method("dist", "optics")
#'
#' d <- dist(random.robinson(50, pre=TRUE, noise=.1))
#'
#' o <- seriate(d, method = "optics")
#' pimage(d, o)
#' }
#'
#' @export
register_optics <- function() {
  check_installed("dbscan")

  .contr <- structure(
    list(eps = NULL,
         minPts = 5),
    help = list(eps = "upper limit of the size of the epsilon neighborhood (see ? optics)" ,
                minPts = "minimum density for dense neighborhoods")
  )

  optics_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    control$minPts <- min(control$minPts, attr(x, "Size"))

    dbscan::optics(x, eps = control$eps, minPts = control$minPts)$order
  }

  set_seriation_method(
    "dist",
    "optics",
    optics_order,
    "Use ordering points to identify the clustering structure (OPTICS) to create an order",
    .contr,
    verbose = TRUE
  )
}
