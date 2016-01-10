#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
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


## registers seriation methods and criteria from package DendSer

register_DendSer <- function() {
  if(!.installed("DendSer")) stop("Package 'DendSer' needs to be  installed!")

  ## seriation methods

  ## control:
  # cost (default: costBAR)
  #   ## costLS, costPL, costLPL, costED, costARc, costBAR
  # h (default is NULL -> complete)

  DendSer_helper <- function(x, control) {
    n <- attr(x, "Size")

    control <- .get_parameters(control, list(
      h = NULL,
      method = "complete",
      criterion = NULL,
      cost = DendSer::costBAR,
      DendSer_args = NULL,
      verbose = FALSE
    ))

    ## fix cost if it is a criterion from seriation
    if(!is.null(control$criterion))
      control$cost <- DendSer::crit2cost(crit = control$criterion)

    ## produce hclust
    if(is.null(control$h))
      control$h <- hclust(x, control$method)

    control$method <- NULL
    control$criterion <- NULL
    control$ser_weight <- x

    if(!is.null(control$DendSer_args)) {
      control <- c(control, control$DendSer_args)
      control$DendSer_args <- NULL
    }

    do.call(DendSer::DendSer, control)
  }


  DendSer_BAR <- DendSer_helper

  DendSer_PL <- function(x, control) {
    control$cost <- DendSer::costPL
    DendSer_helper(x, control)
  }

  DendSer_ARc <- function(x, control) {
    control$cost <- DendSer::costARc
    DendSer_helper(x, control)
  }

  DendSer_LS <- function(x, control) {
    control$cost <- DendSer::costLS
    control$h <- hclust(x)
    DendSer_helper(as.matrix(x)[,1], control)
  }

  seriation::set_seriation_method("dist", "DendSer",
    DendSer_helper, "Dendrogram seriation (DendSer)")

  seriation::set_seriation_method("dist", "DendSer_BAR",
    DendSer_BAR, "Dendrogram seriation (BAR)")
  seriation::set_seriation_method("dist", "DendSer_PL",
    DendSer_PL, "Dendrogram seriation (Path length)")
  seriation::set_seriation_method("dist", "DendSer_ARc",
    DendSer_ARc, "Dendrogram seriation (ARc)")
  seriation::set_seriation_method("dist", "DendSer_LS",
    DendSer_LS, "Dendrogram seriation (Leaf sort)")


  ## criteria

  DendSer_crit_ARc <- function(x, order, ...) {
    x <- as.matrix(x)
    if (is.null(order)) order <- 1:nrow(x)
    else order <- get_order(order)
    DendSer::costARc(x,order,...)
  }

  seriation::set_criterion_method("dist", "ARc", DendSer_crit_ARc,
    "AR cost", FALSE)

  ## Already in seriation
  #  DendSer_crit_BAR <- function(x, order, ...) {
  #    x <- as.matrix(x)
  #    if (is.null(order)) order <- 1:nrow(x)
  #    else order <- get_order(order)
  #    DendSer::costBAR(x,order,...)
  #  }
  #
  #  seriation::set_criterion_method("dist", "BAR", DendSer_crit_BAR,
  #    "Banded AR cost", FALSE)


  #  criterion_method_dist_LPL <- function(x, order, ...) {
  #    x <- as.matrix(x)
  #    if (is.null(order)) order <- 1:nrow(x)
  #    else order <- get_order(order)
  #    DendSer::costLPL(x,order,...)
  #  }
  #
  #  seriation::set_criterion_method("dist", "LPL", criterion_method_dist_LPL,
  #    "Lazy path cost", FALSE)
  #}
}
