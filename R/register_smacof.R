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
#' Registers the `"MDS_smacof"` method for [seriate()] based on multidemensional
#' scaling using stress majorization and the corresponding `"smacof_stress0"`
#' criterion implemented in package smacof (de Leeuw & Mair, 2009).
#'
#' Seriation method `"smacof"` implements stress majorization with several transformation functions.
#' These functions are passed on as the type control parameter. We default
#' to `"ratio"`, which together with `"interval"` performs metric MDS.
#' `"ordinal"` can be used
#' for non-metric MDS. See [smacof::smacofSym()] for details on the
#' control parameters.
#'
#' The corresponding criterion calles `"smacof_stress0"` is also registered.
#' There additional parameter `type` is used to specify the used
#' transformation function. It should agree with the function used for seriation.
#' See [smacof::stress0()] for details on the stress calculation.
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
#' get_seriation_method("dist", "MDS_smacof")
#'
#' d <- dist(random.robinson(20, pre = TRUE))
#'
#' ## use Banded AR form with default clustering (complete-link)
#' o <- seriate(d, "MDS_smacof", verbose = TRUE)
#' pimage(d, o)
#'
#' # recalculate stress for the order
#' MDS_stress(d, o)
#'
#' # ordinal MDS. stress needs to be calculated using the correct type with stress0
#' o <- seriate(d, "MDS_smacof", type = "ordinal", verbose = TRUE)
#' criterion(d, o, method = "smacof_stress0", type = "ordinal")
#' }
#' @export
register_smacof <- function() {
  check_installed("smacof")

  .smacof_control <- structure(
    list(
      type = "ratio",
      init = "torgerson",
      relax = FALSE,
      modulus = 1,
      itmax = 1000,
      eps = 1e-06,
      verbose = FALSE
    ),
    help = list(
      type = 'MDS type: "interval", "ratio", "ordinal" (nonmetric MDS)',
      init = 'start configuration method ("torgerson"/"random")',
      relax = "use block relaxation for majorization?",
      modulus = "number of smacof iterations per monotone regression call",
      itmax = "maximum number of iterations",
      eps = "convergence criterion"
    )
  )

  seriate_dist_smacof <- function(x, control = NULL) {
    control <- .get_parameters(control, .smacof_control)

    r <-
      smacof::smacofSym(
        x,
        ndim = 1,
        type = control$type,
        verbose = control$verbose,
        init = control$init,
        relax = control$relax,
        modulus = control$modulus,
        itmax = control$itmax,
        eps = control$eps
      )

    if (control$verbose)
      print(r)

    config <- drop(r$conf)
    names(config) <- labels(x)

    o <- order(config)

    attr(o, "configuration") <- config
    o
  }


  seriation::set_seriation_method(
    "dist",
    "MDS_smacof",
    seriate_dist_smacof,
    "Seriation based on multidemensional scaling using stress majorization (de Leeuw & Mair, 2009).",
    .smacof_control,
    optimizes = "Other (MDS stress)",
    verbose = TRUE
  )

  smacof_crit_stress0 <-
    function(x,
             order,
             type = "ratio",
             warn = FALSE,
             ...) {
      conf <- get_config(order)
      if (is.null(conf))
        conf <- uniscale(x, order, warn = warn)

      smacof::stress0(x, cbind(conf), type = type, ...)$stress
    }

  seriation::set_criterion_method(
    "dist",
    "smacof_stress0",
    smacof_crit_stress0,
    "Stress0 calculated for different transformation types from package smacof.",
    FALSE,
    verbose = TRUE
  )

}
