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
#' Registers the method `"tsne"` for [seriate()]. This method applies
#' 1D t-SNE to a data matrix or a distance matrix and extracts the order
#' from the 1D embedding. To speed up the process, an initial embedding is
#' created using 1D multi-dimensional scaling (MDS) or principal
#' comonents analysis (PCA) which is improved by t-SNE.
#'
#' The `control` parameter `"mds"` or `"pca"` controls if MDS (for distances)
#' or PCA (for data matrices) is used to create an
#' initial embedding. See [Rtsne::Rtsne()] to learn about the other
#' available `control` parameters.
#'
#' Perplexity is automatically set as the minimum between 30 and the number of
#' observations. It can be also specified using the control parameter
#' `"preplexity"`.
#'
#' **Note:** Package \pkg{Rtsne} needs to be installed.
#'
#' @aliases register_tsne tsne tSNE
#' @seealso [Rtsne::Rtsne()]
#' @family seriation
#' @returns Nothing.
#'
#' @references van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing
#' High-Dimensional Data Using t-SNE. _Journal of Machine Learning Research,_
#' **9**,
#' pp.2579-2605.
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_tsne()
#'
#' # distances
#' get_seriation_method("dist", "tsne")
#'
#' data(SupremeCourt)
#' d <- as.dist(SupremeCourt)
#'
#' o <- seriate(d, method = "tsne", verbose = TRUE)
#' pimage(d, o)
#'
#' # look at the returned configuration and plot it
#' attr(o[[1]], "configuration")
#' plot_config(o)
#'
#' # the t-SNE results are also available as an attribute (see ? Rtsne::Rtsne)
#' attr(o[[1]], "model")
#'
#' ## matrix
#' get_seriation_method("matrix", "tsne")
#'
#' data("Zoo")
#' x <- Zoo
#'
#' x[,"legs"] <- (x[,"legs"] > 0)
#'
#' # t-SNE does not allow duplicates
#' x <- x[!duplicated(x), , drop = FALSE]
#'
#' class <- x$class
#' label <- rownames(x)
#' x <- as.matrix(x[,-17])
#'
#' o <- seriate(x, method = "tsne", eta = 10, verbose = TRUE)
#' pimage(x, o, prop = FALSE, row_labels = TRUE, col_labels = TRUE)
#'
#' # look at the row embedding
#' plot_config(o[[1]], col = class)
#' }
#'
#' @export
register_tsne <- function() {
  check_installed("Rtsne")

  .contr <- structure(
    list(
      max_iter = 1000,
      theta = 0.5,
      perplexity = NULL,
      eta = 100,
      mds = TRUE,
      verbose = FALSE
    ),
    help = list(
      max_iter = "number of iterations",
      theta = "speed/accuracy trade-off (increase for less accuracy)",
      perplexity = "perplexity parameter (calculated as n - 1 / 3)",
      eta = "learning rate",
      mds = "start from a classical MDS solution"
    )
  )

  tsne_order <- function(x, control) {
    control <- .get_parameters(control, .contr)

    # start with MDS
    if (control$mds)
      Y_init <- stats::cmdscale(x, k = 1)
    else
      Y_init <- NULL

    # default is 30 (reduced for low n)
    if (is.null(control$preplexity))
      control$perplexity <- 30

    control$perplexity <-
      max(min(control$perplexity, floor(attr(x, "Size") / 3) - 1), 1)

    embedding <- Rtsne::Rtsne(
      x,
      dims = 1,
      is_distance = TRUE,
      max_iter = control$max_iter,
      theta = control$theta,
      eta = control$eta,
      perplexity = control$perplexity,
      Y_init = Y_init,
      verbose = control$verbose
    )

    o <- order(embedding$Y)

    attr(o, "configuration") <-
      structure(drop(embedding$Y), names = attr(x, "Labels"))
    attr(o, "model") <-  embedding

    o
  }

  .contr_matrix <- structure(
    list(
      max_iter = 1000,
      theta = 0.5,
      perplexity = NULL,
      eta = 100,
      pca = TRUE
    ),
    help = list(max_iter = "number of iterations",
    theta = "speed/accuracy trade-off (increase for less accuracy)",
    perplexity = "perplexity parameter (calculated as n - 1 / 3)",
    eta = "learning rate",
    pca = "start the PCA solution"
  ))

tsne_order_matrix <- function(x, control) {
  control <- .get_parameters(control, .contr_matrix)

  # default is 30 (reduced for low n)
  if (is.null(control$preplexity))
    control$perplexity <- 30

  control$perplexity <-
    max(min(control$perplexity, floor(nrow(x) / 3) - 1), 1)

  embedding <-
    Rtsne::Rtsne(
      x,
      dims = 1,
      is_distance = FALSE,
      pca = control$pca,
      max_iter = control$max_iter,
      theta = control$theta,
      eta = control$eta,
      perplexity = control$perplexity,
      verbose = control$verbose
    )

  o <- order(embedding$Y)

  attr(o, "configuration") <-
    structure(drop(embedding$Y), names = rownames(x))
  attr(o, "model") <-  embedding

  o
}

tsne_order_matrix_2 <-
  function(x, control, margin = seq_along(dim(x))) {
    if (1L %in% margin)
      row <- tsne_order_matrix(x, control)
    else
      row <- NA

    if (2L %in% margin)
      col <- tsne_order_matrix(t(x), control)
    else
      col <- NA

    list(row, col)
  }

set_seriation_method(
  "dist",
  "tsne",
  tsne_order,
  "Use 1D t-distributed stochastic neighbor embedding (t-SNE) a distance matrix to create an order (van der Maaten and Hinton, 2008).",
  .contr,
  randomized = TRUE,
  verbose = TRUE
)

set_seriation_method(
  "matrix",
  "tsne",
  tsne_order_matrix_2,
  "Use 1D t-distributed stochastic neighbor embedding (t-SNE) of the rows of a matrix to create an order (van der Maaten and Hinton, 2008).",
  .contr_matrix,
  randomized = TRUE,
  verbose = TRUE
)
}
