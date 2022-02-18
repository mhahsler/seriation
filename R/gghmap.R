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

#' @rdname hmap
#' @include hmap.R
#' @export
gghmap <- function(x,
  distfun = stats::dist,
  method = "OLO",
  control = NULL,
  scale = c("none", "row", "column"),
  prop = FALSE,
  ...) {
  scale <- match.arg(scale)

  if (inherits(x, "dist")) {
    # scale and distFun are ignored!
    o <- seriate(x, method = method, control = control)
  } else {
    x <- as.matrix(x)

    contr <- list(
      dist_fun = distfun,
      seriation_method = method,
      seriation_control = control,
      scale = scale
    )

    o <-
      seriate(x,
        method = "heatmap",
        control = contr)
  }

  ggpimage(x, o, prop = prop, ...)

}
