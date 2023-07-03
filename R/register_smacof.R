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


#' Register Seriation Methods from Package smacof
#'
#' Registers methods for [seriate()] based on multidemensional scaling using
#' stress majorization implemented in package smacof (de Leeuw & Mair, 2009).
#'
#' **Rewite**
#' Registers the method \code{"DendSer"} for [seriate()]. DendSer is a fast
#' heuristic for reordering dendrograms developed by Earle and Hurley (2015)
#' able to use different criteria.
#'
#' \code{control} for \code{seriate} with
#' method \code{"DendSer"} accepts the following parameters:
#'
#' - "h" or "method" A dendrogram or a method for hierarchical clustering
#'   (see \code{hclust}). Default: complete-link.
#' - "criterion" A seriation criterion to optimize (see
#'   \code{list_criterion_methods("dist")}). Default: \code{"BAR"} (Banded
#'   anti-Robinson from with 20\% band width).}
#' - "verbose" a logical; print progress information?
#' - "DendSer_args" additional arguments for \code{DendSer}.
#'
#'
#' Note: Package \pkg{smacof} needs to be installed.
#'
#' @aliases registersmacof smacof
#' @family seriation
#' @returns Nothing.
#'
#' @references
#' Jan de Leeuw, Patrick Mair (2009). Multidimensional Scaling Using Majorization: SMACOF in R.
#' _Journal of Statistical Software, 31(3),_ 1-30. \doi{10.18637/jss.v031.i03}
#' @keywords optimize cluster
#' @examples
#' \dontrun{
#' register_smacof()
#'
#' get_seriation_method("dist", "nMDS")
#'
#' d <- dist(random.robinson(20, pre = TRUE))
#'
#' ## use Banded AR form with default clustering (complete-link)
#' o <- seriate(d, "nMDS", verbose = TRUE)
#' pimage(d, o)
#'
#' # recalculate stress for the order
#' MDS_stress(d, o)
#' }
#' @export
register_smacof <- function() {
  check_installed("smacof")

  .smacof_control <- list(
    init = "torgerson",
    type = "interval",
    relax = FALSE,
    modulus = 1,
    itmax = 1000,
    eps = 1e-06,
    verbose = FALSE
    )

  seriate_dist_smacof <- function(x, control = NULL) {
    control <- .get_parameters(control, .smacof_control)

    # this is nMDS
    r <- smacof::smacofSym(x, ndim = 1, type = control$type, verbose = control$verbose,
                           init = control$init, relax = control$relax,
                           modulus = control$modulus, itmax = control$itmax,
                           eps = control$eps)

    if (control$verbose)
      print(r)

    config <- drop(r$conf)
    names(config) <- labels(x)

    o <- order(config)

    attr(o, "configuration")
    o
  }


  seriation::set_seriation_method(
    "dist",
    "MDS_smacof",
    seriate_dist_smacof,
    "Seriation based on multidemensional scaling using stress majorization using smacof::smacofSym() (de Leeuw & Mair, 2009).",
    .smacof_control,
    optimizes = "Other (MDS stress)"
  )

}
