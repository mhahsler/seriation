#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2017 Michael Hahsler, Christian Buchta and Kurt Hornik
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

# unidimensional scaling: Defrays Decomposition (1978)

#' Unidimensional Scaling from Seriation Results
#'
#' Calculates an (approximate) unidimensional scaling configuration given an order.
#'
#' This implementation uses the method describes in Maier and De Leeuw (2015) to calculate the
#' minimum stress configuration for a given (seriation) order by performing a 1D MDS fit.
#' If the 1D MDS fit does not preserve the given order perfectly (i.e., because the
#' order was not calculated using MDS), then a warning is produced indicating
#' for how many positions order could not be preserved.
#'
#' If no order is specified, then
#' a seriation order is first found using the specified seriation method.
#'
#' The code is similar to `uniscale()` in \pkg{smacof} (de Leeuw, 2090), but scales to larger
#' datasets since it only checks the permutation given by the seriation order.
#'
#' @param d a dissimilarity matrix.
#' @param order a precomputed permutation (configuration) order.  If
#' \code{NULL}, then seriation is performed using the method specified in
#' \code{method}.
#' @param method seriation method used if \code{o} is \code{NULL}.
#' @param rep Number of repetitions of the seriation heuristic.
#' @param warn logical; produce a warning if the 1D MDS fit does not preserve the
#'  given order.
#' @param \dots additional arguments are passed on to the seriation method.
#' @return A vector with the fitted configuration.
#' @author Michael Hahsler with code from Patrick Mair (from \pkg{smacof}).
#' @references Mair P., De Leeuw J. (2015). Unidimensional scaling. In
#' \emph{Wiley StatsRef: Statistics Reference Online,} Wiley, New York.
#' \doi{10.1002/9781118445112.stat06462.pub2}
#'
#' Jan de Leeuw, Patrick Mair (2009). Multidimensional Scaling Using Majorization:
#' SMACOF in R. Journal of Statistical Software, 31(3), 1-30.
#' \doi{10.18637/jss.v031.i03}
#'
#' @keywords optimize
#' @examples
#' data(SupremeCourt)
#' d <- as.dist(SupremeCourt)
#'
#' # embedding-based methods return a configuration as the attribute "embedding"
#' # configplot visualizes the embedding
#' o <- seriate(d, method = "MDS_sammon")
#' get_order(o)
#' attr(o[[1]], "embedding")
#' configplot(o)
#'
#' # angle methods return a 2D configuration
#' o <- seriate(d, method = "MDS_angle")
#' get_order(o)
#' attr(o[[1]], "embedding")
#' configplot(o)
#'
#' # calculate a configuration for a seriation method that does not
#' # use an embedding
#' o <- seriate(d, method = "spectral")
#' get_order(o)
#' attr(o[[1]], "embedding")
#'
#' # find the minimum-stress configuration
#' sc <- uniscale(d, o)
#' sc
#'
#' configplot(sc)
#'
#' @export
uniscale <-
  function(d,
    order = NULL,
    method = "MDS",
    rep = 10,
    warn = TRUE,
    ...) {
    if (is.null(order))
      order <- seriate(d, method = method, rep = rep, ...)

    o <- get_rank(order)
    n <- length(o)

    # we do not use weights
    w <- 1 - diag(n)

    normDissN <- function (diss, wghts, m) {
      N <- length(diss) * m
      dissnorm <- diss / sqrt(sum(wghts * diss ^ 2, na.rm = TRUE)) *
        sqrt(N)
      return(dissnorm)
    }

    delta <- as.matrix(normDissN(as.dist(d), as.dist(w),
      1))
    v <- as.matrix(solve((diag(rowSums(
      w
    )) - w) + (1 / n)) - (1 / n))
    s <- sign(outer(o, o, "-"))

    t <- as.vector(v %*% rowSums(delta * w * s))
    names(t) <- names(o)

    # does the configuration preserve the order in o?
    mismatches <- sum(order(o) != order(t))
    if (mismatches > 0 && warn)
      warning("Configutation order does not preserve given order! Mismatches: ", mismatches, " of ", n)

    t
  }

#' @rdname uniscale
#' @param x a scaling returned by `uniscale()` or a
#'   `ser_permutation` with an embedding attribute.
#' @param main main plot label
#' @param pch print character
#' @export
configplot <- function (x, main, pch = 19, ...) {
  if (missing(main))
    main <- "Configuration"

  if(inherits(x, "ser_permutation"))
    x <- x[[1]]

  if (inherits(x, "ser_permutation_vector")) {
    o <- get_order(x)  # only used for 2D case
    if(!is.null(attr(x, "configuration")))
      x <- attr(x, "configuration")
    else if(!is.null(attr(x, "embedding")))
      x <- attr(x, "embedding")
    else
      stop("Permutation vector has no configuration attribute. Use uniscale() first to calcualte a configuration")
  }

  # 2D
  if (is.matrix(x)) {
    graphics::plot(x, pch = pch, main = main, ...)
    graphics::text(x = x, labels = rownames(x), pos = 1)
    graphics::lines(x[o, , drop = FALSE], col = "grey")
    return()
  }

  # 1D
  x <- drop(x)
  n <- length(x)
  plot(
    x,
    rep(0, n),
    axes = FALSE,
    ann = FALSE,
    pch = pch,
    type = "o",
    ylim = c(-0.2, 0.8),
    ...
  )
  title(main)

  labs <- names(x)
  if (is.null(labs))
    labs <- 1:n
  text(x,
    rep(0, n) + 0.05,
    labs,
    srt = 90,
    adj = c(0, 0.5))
}
