#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2017 Michael Hahsler, Christian Buchta and Kurt Hornik
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


orderplot <- function (x, main, pch = 19, ...){
  if (missing(main)) main <- "Configuration"

  n <- length(x)
  plot(x, rep(0, n), axes = FALSE, ann = FALSE,
    pch = pch, type = "o", ylim = c(-0.2, 0.8), ...)
  title(main)

  labs <- names(x)
  if(is.null(labs)) labs <- 1:n
  text(x, rep(0, n) + 0.05, labs,
    srt = 90, adj = c(0, 0.5))
}

uniscale <- function(d, order = NULL, method = "QAP_LS", rep = 10, ...) {
  if(is.null(order)) order <- seriate(d, method = method, rep = rep, ...)
  x <- get_rank(order)
  n <- length(x)

  w <- 1 - diag(n)

  normDissN <- function (diss, wghts, m){
    N <- length(diss) * m
    dissnorm <- diss/sqrt(sum(wghts * diss^2, na.rm = TRUE)) *
      sqrt(N)
    return(dissnorm)
  }

  delta <- as.matrix(normDissN(as.dist(d), as.dist(w),
    1))
  v <- as.matrix(solve((diag(rowSums(w)) - w) + (1/n)) - (1/n))
  s <- sign(outer(x, x, "-"))

  t <- as.vector(v %*% rowSums(delta * w * s))

  names(t) <- attr(d, "Labels")
  t
}
