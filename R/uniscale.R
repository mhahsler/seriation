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
#' Performs (approximate) unidimensional scaling by first performing seriation
#' to obtain a permutation and the using the permutation to calculate the
#' configuration.
#'
#' Uses the method describes in Maier and De Leeuw (2015) to calculate the
#' minimum stress configuration for either a given
#' configuration/permutation/order or for a permutation computed via a
#' seriation method.
#'
#' The code is similar to `uniscale()` in \pkg{smacof}, but scales to larger
#' datasets since it does not check all permutations.
#'
#' @param d a dissimilarity matrix.
#' @param order a precomputed permutation (configuration) order.  If
#' \code{NULL}, then seriation is performed using the method specified in
#' \code{method}.
#' @param method seriation method used if \code{o} is \code{NULL}.
#' @param rep Number of repetitions of the seriation heuristic.
#' @param \dots additional arguments are passed on to the seriation method.
#' @return A vector with the fitted configuration.
#' @author Michael Hahsler with code from Patrick Mair (from \pkg{smacof}).
#' @references Mair P., De Leeuw J. (2015). Unidimensional scaling. In
#' \emph{Wiley StatsRef: Statistics Reference Online,} Wiley, New York.
#' \doi{10.1002/9781118445112.stat06462.pub2}
#' @keywords optimize
#' @examples
#' data(SupremeCourt)
#'
#' d <- as.dist(SupremeCourt)
#'
#' sc <- uniscale(d)
#' sc
#'
#' orderplot(sc)
#' @export
uniscale <-
  function(d,
    order = NULL,
    method = "QAP_LS",
    rep = 10,
    ...) {
    if (is.null(order))
      order <- seriate(d, method = method, rep = rep, ...)
    x <- get_rank(order)
    n <- length(x)

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
    s <- sign(outer(x, x, "-"))

    t <- as.vector(v %*% rowSums(delta * w * s))

    names(t) <- attr(d, "Labels")
    t
  }

#' @rdname uniscale
#' @param x a scaling returned by `uniscale()`.
#' @param main main plot label
#' @param pch print character
#' @export
orderplot <- function (x, main, pch = 19, ...) {
  if (missing(main))
    main <- "Configuration"

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
