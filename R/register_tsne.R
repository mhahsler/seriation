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

#' Register Seriation Based on 1D t-SNE
#'
#' Use t-distributed stochastic neighbor embedding (t-SNE) for [seriate()].
#'
#' Registers the method \code{"tsne"} for [seriate()]. This method applies
#' 1D t-SNE to data represented by a distance matrix and extracts the order
#' from the 1D embedding. To speed up the process, an initial embedding is
#' created using multi-dimensional scaling (MDS) which is improved by t-SNE.
#'
#' The \code{control} parameter \code{mds} controls if MDS is used to create an
#' initial embedding. See [Rtsne::Rtsne()] to learn about the other
#' available \code{control} parameters.
#'
#' \bold{Note:} Package \pkg{Rtsne} needs to be installed.
#'
#' @aliases register_tsne tsne tSNE
#' @seealso [Rtsne::Rtsne()] in \pkg{Rtsne}.
#' @family seriation
#' @returns Nothing.
#'
#' @references van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing
#' High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 9,
#' pp.2579-2605.
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_tsne()
#' get_seriation_method("dist", "tsne")
#'
#' d <- dist(random.robinson(50, pre=TRUE, noise=.1))
#'
#' o <- seriate(d, method = "tsne")
#' pimage(d, o)
#' }
#'
register_tsne <- function() {
  check_installed("Rtsne")

  .contr <- list(
    max_iter = 1000,
    theta = 0,
    perplexity = 30,
    eta = 200,
    mds = TRUE
  )

  tsne_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    # start with MDS
    if(control$mds) Y_init <- stats::cmdscale(x, k = 1)
    else Y_init <- NULL

    # calculate the maximal value for perplexity
    perplexity <- min(control$perplexity, floor(attr(x, "Size") / 3) - 1)

    embedding <- Rtsne::Rtsne(x, dims = 1, is_distance = TRUE,
      max_iter = control$max_iter, theta = control$theta, eta = control$eta,
      perplexity = perplexity, Y_init = Y_init)
    order(embedding$Y)
  }

  set_seriation_method(
    "dist",
    "tsne",
    tsne_order,
    "Use 1D  t-distributed stochastic neighbor embedding (t-SNE) to one dimension to create an order",
    .contr
  )
}
