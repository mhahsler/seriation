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

## register seriation based on t-SNE

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
    if(control$mds) Y_init <- cmdscale(x, k = 1)
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
