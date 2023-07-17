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
#' Registers the method `"umap"` for [seriate()]. This method applies
#' 1D UMAP to a data matrix or a distance matrix and extracts the order from
#' the 1D embedding.
#'
#' Control parameter `n_epochs` can be increased to find a better embedding.
#'
#' The returned seriation permutation vector has an attribute named
#' `embedding` containing the umap embedding.
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
#'
#' ## distances
#' get_seriation_method("dist", "umap")
#'
#' data(SupremeCourt)
#' d <- as.dist(SupremeCourt)
#'
#' o <- seriate(d, method = "umap", verbose = TRUE)
#' pimage(d, o)
#'
#' # look at the returned embedding and plot it
#' attr(o[[1]], "configuration")
#' plot_config(o)
#'
#' ## matrix
#' get_seriation_method("matrix", "umap")
#'
#' data("Zoo")
#' Zoo[,"legs"] <- (Zoo[,"legs"] > 0)
#' x <- as.matrix(Zoo[,-17])
#' label <- rownames(Zoo)
#' class <- Zoo$class
#'
#' o <- seriate(x, method = "umap", verbose = TRUE)
#' pimage(x, o)
#'
#' plot_config(o[[1]], col = class)
#' }
#' @export
register_umap <- function() {
  check_installed("umap")

  .contr <- unclass(umap::umap.defaults)
  .contr$n_epochs <- 1000
  .contr$n_neighbors <- NA
  .contr$n_components <- 1
  .contr$alpha <- 0.001
  .contr$input <- NA
  .contr$random_state <- NA

  attr(.contr, "help") <- list(n_neighbors = "see ? umap::umap for help")


  umap_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    if (is.na(control$input))
      control$input <- "dist"

    x <- as.matrix(x)

    # we cannot have more neighbors than data points
    if (is.na(control$n_neighbors))
      control$n_neighbors <- 15
    control$n_neighbors <- min(control$n_neighbors, nrow(x))

    # use different random numbers for every run
    if (is.na(control$random_state))
      control$random_state <-
      as.integer(runif(1, 0, .Machine$integer.max))

    # has to be 1
    control$n_components <- 1

    class(control) <- class(umap::umap.defaults)

    embedding <- umap::umap(x, config = control)

    o <- order(embedding$layout)
    embedding <- drop(embedding$layout)
    names(embedding) <- rownames(x)
    attr(o, "configuration") <- embedding

    o
  }

  umap_order_matrix_2 <-
    function(x, control, margin = seq_along(dim(x))) {
      control$input <- "data"

      if (1L %in% margin)
        row <- umap_order(x, control)
      else
        row <- NA

      if (2L %in% margin)
        col <- umap_order(t(x), control)
      else
        col <- NA

      list(row, col)
    }

  set_seriation_method(
    "dist",
    "umap",
    umap_order,
    "Use 1D Uniform manifold approximation and projection (UMAP) embedding of the distances to create an order",
    .contr,
    randomized = TRUE,
    verbose = TRUE
  )

  set_seriation_method(
    "matrix",
    "umap",
    umap_order_matrix_2,
    "Use 1D Uniform manifold approximation and projection (UMAP) embedding of the data to create an order",
    .contr,
    randomized = TRUE,
    verbose = TRUE
  )
}
