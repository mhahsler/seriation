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


#' Register Seriation Methods from Package DendSer
#'
#' Register the DendSer dendrogram seriation method and the ARc criterion
#' (Earle and Hurley, 2015) for use with [seriate()].
#'
#' Registers the method `"DendSer"` for seriate. DendSer is a fast
#' heuristic for reordering dendrograms developed by Earle and Hurley (2015)
#' able to use different criteria.
#'
#' `control` for [`seriate()`] with
#' method `"DendSer"` accepts the following parameters:
#'
#' - `"h"` or `"method"`: A dendrogram or a method for hierarchical clustering
#'   (see [hclust]). Default: complete-link.
#' - `"criterion"`: A seriation criterion to optimize (see
#'   `list_criterion_methods("dist")`. Default: `"BAR"` (Banded
#'   anti-Robinson from with 20% band width).
#' - `"verbose"`: a logical; print progress information?
#' - `"DendSer_args"`: additional arguments for [`DendSer::DendSer()`].
#'
#' For convenience, the following methods (for different cost functions) are
#' also provided:
#'
#' - `"DendSer_ARc"` (anti-robinson form),
#' - `"DendSer_BAR"` (banded anti-Robinson form),
#' - `"DendSer_LPL"` (lazy path length),
#' - `"DendSer_PL"` (path length).
#'
#' **Note:** Package \pkg{DendSer} needs to be installed.
#'
#' @aliases register_DendSer DendSer dendser
#' @seealso [`DendSer::DendSer()`]
#' @family seriation
#' @returns Nothing.
#'
#' @author Michael Hahsler based on code by Catherine B. Hurley and Denise
#' Earle
#' @references D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation
#' for application to visualization. _Journal of Computational and
#' Graphical Statistics,_ **24**(1), 1--25.
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_DendSer()
#' get_seriation_method("dist", "DendSer")
#'
#' d <- dist(random.robinson(20, pre=TRUE))
#'
#' ## use Banded AR form with default clustering (complete-link)
#' o <- seriate(d, "DendSer_BAR")
#' pimage(d, o)
#'
#' ## use different hclust method (Ward) and AR as the cost function for
#' ## dendrogram reordering
#' o <- seriate(d, "DendSer", control = list(method = "ward.D2", criterion = "AR"))
#' pimage(d, o)
#' }
#'
#' @export
register_DendSer <- function() {
  check_installed("DendSer")

  ## seriation methods

  ## control:
  # cost (default: costBAR)
  #   ## costLS, costPL, costLPL, costED, costARc, costBAR
  # h (default is NULL -> complete)


  .DendSer_control <- structure(
    list(
      h = NULL,
      method = "complete",
      criterion = NULL,
      DendSer_args = NULL,
      verbose = FALSE
    ),
    help = list(
      h = "an hclust object (optional)",
      method = "hclust linkage method",
      criterion = "criterion to optimize the dendrogram for",
      DendSer_args = "more arguments are passed on to DendSer (? DendSer)"
    )

  )

  DendSer_helper <- function(x, control) {
    n <- attr(x, "Size")

    control <- .get_parameters(control, .DendSer_control)

    control$cost <- DendSer::crit2cost(crit = control$criterion)
    control$criterion <- NULL

    ## produce hclust
    if (is.null(control$h))
      control$h <- hclust(x, control$method)
    control$method <- NULL

    control$ser_weight <- x

    if (!is.null(control$DendSer_args)) {
      control <- c(control, control$DendSer_args)
      control$DendSer_args <- NULL
    }

    permute(control$h, do.call(DendSer::DendSer, control))
  }


  DendSer_BAR <-  function(x, control) {
    control$criterion <- "BAR"
    DendSer_helper(x, control)
  }


  DendSer_PL <- function(x, control) {
    control$criterion <- "Path_length"
    DendSer_helper(x, control)
  }

  DendSer_LPL <- function(x, control) {
    control$criterion <- "Lazy_path_length"
    DendSer_helper(x, control)
  }

  DendSer_ARc <- function(x, control) {
    control$criterion <- "Arc"
    DendSer_helper(x, control)
  }

  ## This is not Least Squares!
  #  DendSer_LS <- function(x, control) {
  #    control$cost <- DendSer::costLS
  #    control$criterion <- "LS"
  #    control$h <- hclust(x)
  #    DendSer_helper(as.matrix(x)[,1], control)
  #  }

  set_seriation_method(
    "dist",
    "DendSer",
    DendSer_BAR,
    "Dendrogram seriation (Earle and Hurley, 2015).",
    .DendSer_control,
    optimizes = .opt(NA, "specified criterion restricted by dendrogram"),
    verbose = TRUE
  )

  set_seriation_method(
    "dist",
    "DendSer_BAR",
    DendSer_BAR,
    "Dendrogram seriation with BAR (Earle and Hurley, 2015).",
    .DendSer_control,
    optimizes = .opt("BAR", "banded anti-Robinson form restricted by dendrogram"),
    verbose = TRUE
  )

  set_seriation_method(
    "dist",
    "DendSer_PL",
    DendSer_PL,
    "Dendrogram seriation for Path length (Earle and Hurley, 2015).",
    .DendSer_control,
    optimizes = .opt("Path_length", "restricted by dendrogram"),
    verbose = TRUE
  )

  set_seriation_method(
    "dist",
    "DendSer_LPL",
    DendSer_LPL,
    "Dendrogram seriation for Lazy path length (Earle and Hurley, 2015).",
    .DendSer_control,
    optimizes = .opt("Lazy_path_length", "restricted by dendrogram"),
    verbose = TRUE
  )

  set_seriation_method(
    "dist",
    "DendSer_ARc",
    DendSer_ARc,
    "Dendrogram seriation for Anti-Robinson form cost (Earle and Hurley, 2015).",
    optimizes = .opt("ARc", "Anti-Robinson form cost restricted by dendrogram"),
    .DendSer_control,
    verbose = TRUE
  )

  #  set_seriation_method("dist", "DendSer_LS",
  #    DendSer_LS, "Dendrogram seriation (Leaf sort)")

  ## criteria

  DendSer_crit_ARc <- function(x, order, ...) {
    x <- as.matrix(x)
    if (is.null(order))
      order <- 1:nrow(x)
    else
      order <- get_order(order)
    DendSer::costARc(x, order, ...)
  }

  set_criterion_method("dist", "ARc", DendSer_crit_ARc,
                                  "Anti-Robinson form cost (Earle and Hurley, 2015).", FALSE, verbose = TRUE)

  ## Already in seriation
  #  DendSer_crit_BAR <- function(x, order, ...) {
  #    x <- as.matrix(x)
  #    if (is.null(order)) order <- 1:nrow(x)
  #    else order <- get_order(order)
  #    DendSer::costBAR(x,order,...)
  #  }
  #
  #  set_criterion_method("dist", "BAR", DendSer_crit_BAR,
  #    "Banded AR cost", FALSE)


  #  criterion_method_dist_LPL <- function(x, order, ...) {
  #    x <- as.matrix(x)
  #    if (is.null(order)) order <- 1:nrow(x)
  #    else order <- get_order(order)
  #    DendSer::costLPL(x,order,...)
  #  }
  #
  #  set_criterion_method("dist", "LPL", criterion_method_dist_LPL,
  #    "Lazy path cost", FALSE)
  #}
}
