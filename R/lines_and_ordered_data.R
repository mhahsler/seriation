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

#' Create Simulated Data for Seriation Evaluation
#'
#' Several functions to create simulated data to evaluate different aspects of
#' seriation algorithms and criterion functions.
#'
#' `create_lines_data()` creates the lines data set used in for [iVAT()] in
#' Havens and Bezdeck (2012).
#'
#' `create_ordered_data()` is a versatile function which creates "orderable"
#' 2D data using Gaussian components along a linear or circular path. The
#' components are equally spaced (`spacing`) along the path. The default
#' spacing of 6 ensures that 2 adjacent components with a standard deviation of
#' one along the direction of the path will barely touch. The standard
#' deviation along the path is set by `sd1`. The standard deviation
#' perpendicular to the path is set by `sd2`. A value larger than zero
#' will result in the data not being perfectly orderable (i.e., the resulting
#' distance matrix will not be a perfect pre-anti-Robinson matrix and contain
#' anti-Robinson violation events after seriation). Note that a circular path
#' always creates anti-Robinson violation since the circle has to be broken at
#' some point to create a linear order. This function was created for this package
#' (Hahsler et al, 2021).
#'
#' @param n number of data points to create.
#' @param k number of Gaussian components.
#' @param size relative size (number of points) of components (length of k).
#' If `NULL` then all components have the same size.
#' @param spacing space between the centers of components. The default of 6
#' means that the components will barely touch at `ds1 = 1` (3 standard
#' deviations for each Gaussian component).
#' @param path Are the components arranged along a `"linear"` or
#' `"circular"` path?
#' @param sd1 variation in the direction along the components.  A value greater
#' than one means the components are mixing.
#' @param sd2 variation perpendicular to the direction along the components.  A
#' value greater than 0 will introduce anti-Robinson violation events.
#' @returns a data.frame with the created data.
#'
#' @author Michael Hahsler
#' @seealso [seriate()], [criterion()], [iVAT()].
#' @references
#' Havens, T.C. and Bezdek, J.C. (2012): An Efficient Formulation
#' of the Improved Visual Assessment of Cluster Tendency (iVAT) Algorithm,
#' \emph{IEEE Transactions on Knowledge and Data Engineering,} \bold{24}(5),
#' 813--822.
#'
#' Michael Hahsler, Christian Buchta and Kurt Hornik (2021). seriation: Infrastructure for
#' Ordering Objects Using Seriation. R package version 1.3.2.
#' \url{https://github.com/mhahsler/seriation}
#' @keywords datasets
#' @examples
#'
#' ## lines data set from Havens and Bezdek (2011)
#' x <- create_lines_data(250)
#' plot(x, xlim = c(-5, 5), ylim = c(-3, 3), cex = .2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "OLO_single"), col = bluered(100, bias = .5), key = TRUE)
#'
#' ## create_ordered_data can produce many types of "orderable" data
#'
#' ## perfect pre-Anti-Robinson matrix (with a single components)
#' x <- create_ordered_data(250, k = 1)
#' plot(x, cex = .2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "MDS"), col = bluered(100, bias=.5), key = TRUE)
#'
#' ## separated components
#' x <- create_ordered_data(250, k = 5)
#' plot(x, cex =.2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "MDS"), col = bluered(100, bias = .5), key = TRUE)
#'
#' ## overlapping components
#' x <- create_ordered_data(250, k = 5, sd1 = 2)
#' plot(x, cex = .2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "MDS"), col = bluered(100, bias = .5), key = TRUE)
#'
#' ## introduce anti-Robinson violations (a non-zero y value)
#' x <- create_ordered_data(250, k = 5, sd1 = 2, sd2 = 5)
#' plot(x, cex = .2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "MDS"), col = bluered(100, bias = .5), key = TRUE)
#'
#' ## circular path (has always violations)
#' x <- create_ordered_data(250, k = 5, path = "circular", sd1 = 2)
#' plot(x, cex = .2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "OLO"), col = bluered(100, bias = .5), key = TRUE)
#'
#' ## circular path (with more violations violations)
#' x <- create_ordered_data(250, k = 5, path = "circular", sd1 = 2, sd2 = 1)
#' plot(x, cex=.2, col = attr(x, "id"))
#' d <- dist(x)
#' pimage(d, seriate(d, "OLO"), col = bluered(100, bias = .5), key = TRUE)
#' @export
create_lines_data <- function(n = 250) {
  n1 <- n / 5 * 2
  n2 <- n / 5
  n3 <- n / 5 * 2

  x1 <-
    data.frame(x = runif(n1, -5, 5), y = rnorm(n1, mean = 2, sd = .1))
  x2 <-
    data.frame(x = runif(n2, -3, 3), y = rnorm(n2, mean = 0, sd = .1))
  x3 <-
    data.frame(x = runif(n3, -5, 5), y = rnorm(n3, mean = -2, sd = .1))
  id <-
    c(rep(1, times = n1), rep(2, times = n2), rep(3, times = n3))


  x <- rbind(x1, x2, x3)
  o <- sample(nrow(x))
  x <- x[o,]
  id <- id[o]

  rownames(x) <- 1:nrow(x)
  attr(x, "id") <- id

  x
}

#' @rdname create_lines_data
#' @export
create_ordered_data <- function(n = 250,
  k = 2,
  size = NULL,
  spacing = 6,
  path = "linear",
  sd1 = 1,
  sd2 = 0) {
  if (k > n)
    stop("k needs to be less than n!")
  path <- match.arg(path, c("linear", "circular"))

  ## size
  if (is.null(size))
    size <- rep(1, k)
  else if (length(size) != k)
    stop("length of size vector and k do not agree!")
  size <- round(size / sum(size) * n)
  size[1] <- n - sum(size[-1])

  ## create data
  ids <- rep(1:k, times = size)

  x <- data.frame(x = rnorm(n, mean = ids * spacing, sd = sd1),
    y = rnorm(n, mean = 0, sd = sd2))

  ## transform
  if (path == "circular") {
    p <- k * spacing
    theta <- x[, 1] / p * 2 * pi
    r <-  p / (2 * pi) + x[, 2]
    x <- cbind(x = r * sin(theta), y = r * cos(theta))
  }

  ## randomize order
  o <- sample(nrow(x))
  x <- x[o , , drop = FALSE]
  ids <- ids[o]
  attr(x, "id") <- ids

  x
}
