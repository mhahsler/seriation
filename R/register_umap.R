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


#' Register Seriation Based on 1D UMAP
#'
#' Use uniform manifold approximation and projection (UMAP) to embed the data
#' on the number line and create a order for [seriate()].
#'
#' Registers the method \code{"umap"} for [seriate()]. This method applies
#' 1D UMAP to data represented by a distance matrix and extracts the order from
#' the 1D embedding.
#'
#' \bold{Note:} Package \pkg{umap} needs to be installed.
#'
#' @aliases register_umap umap
#' @seealso [umap::umap()] in \pkg{umap}.
#' @family seriation
#' @returns Nothing.
#'
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and
#' Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_umap()
#' get_seriation_method("dist", "umap")
#'
#' d <- dist(random.robinson(50, pre=TRUE, noise=.1))
#'
#' o <- seriate(d, method = "umap")
#' pimage(d, o)
#' }
#'
#' @export
register_umap <- function() {
  check_installed("umap")

  .contr <- unclass(umap::umap.defaults)
  .contr$n_components <- 1
  .contr$input <- "dist"

  umap_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    class(control) <- class(umap::umap.defaults)

    embedding <- umap::umap(as.matrix(x), config = control)
    order(embedding$layout)
  }

  set_seriation_method(
    "dist",
    "umap",
    umap_order,
    "Use 1D Uniform manifold approximation and projection (UMAP) embedding to create an order",
    .contr
  )
}
