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

.dist_methods <-
  c("spearman",
    "kendall",
    "manhattan",
    "euclidean",
    "hamming",
    "ppc",
    "aprd")

#' Dissimilarities and Correlations Between Seriation Orders
#'
#' Calculates dissimilarities/correlations between seriation orders in a list of type
#' [ser_permutation_vector].
#'
#' `ser_cor()` calculates the correlation between two sequences (orders).
#' Note that a seriation order and its reverse are identical and purely an
#' artifact due to the method that creates the order. This is a major
#' difference to rankings. For ranking-based correlation measures (Spearman and
#' Kendall) the absolute value of the correlation is returned for
#' `reverse = TRUE` (in effect returning the correlation for the reversed order). If
#' `test = TRUE` then the appropriate test for association is performed
#' and a matrix with p-values is returned as the attribute `"p-value"`.
#' Note that no correction for multiple testing is performed.
#'
#' For `ser_dist()`, the correlation coefficients (Kendall's tau and
#' Spearman's rho) are converted into a dissimilarity by taking one minus the
#' correlation value. Note that Manhattan distance between the ranks in a
#' linear order is equivalent to Spearman's footrule metric (Diaconis 1988).
#' `reverse = TRUE` returns the pairwise minima using also reversed
#' orders.
#'
#' The positional proximity coefficient (ppc) is a precedence invariant measure
#' based on product of the squared positional distances in two permutations
#' defined as (see Goulermas et al 2016):
#'
#' \deqn{d_{ppc}(R, S) = 1/h \sum_{j=2}^n \sum_{i=1}^{j-1}
#' (\pi_R(i)-\pi_R(j))^2 * (\pi_S(i)-\pi_S(j))^2,}
#'
#' where \eqn{R} and \eqn{S} are two seriation orders, \eqn{pi_R} and
#' \eqn{pi_S} are the associated permutation vectors and \eqn{h} is a
#' normalization factor. The associated generalized correlation coefficient is
#' defined as \eqn{1-d_{ppc}}. For this precedence invariant measure
#' `reverse` is ignored.
#'
#' The absolute pairwise rank difference (aprd) is also precedence invariant
#' and defined as a distance measure:
#'
#' \deqn{d_{aprd}(R, S) = \sum_{j=2}^n \sum_{i=1}^{j-1} | |\pi_R(i)-\pi_R(j)| -
#' |\pi_S(i)-\pi_S(j)| |^p,}
#'
#' where \eqn{p} is the power which can be passed on as parameter `p` and
#' is by default set to 2. For this precedence invariant measure `reverse`
#' is ignored.
#'
#' `ser_align()` tries to normalize the direction in a list of seriations
#' such that ranking-based methods can be used. We add for each permutation
#' also the reversed order to the set and then use a modified version of Prim's
#' algorithm for finding a minimum spanning tree (MST) to choose if the
#' original seriation order or its reverse should be used. We use the orders
#' first added to the MST. Every time an order is added, its reverse is removed
#' from the possible remaining orders.
#'
#' @family permutation
#'
#' @param x set of seriation orders as a list with elements which can be
#' coerced into [ser_permutation_vector] objects.
#' @param y if not `NULL` then a single seriation order can be specified.
#' In this case `x` has to be a single seriation order and not a list.
#' @param method a character string with the name of the used measure.
#' Available measures are: `"kendall"`, `"spearman"`,
#' `"manhattan"`, `"euclidean"`, `"hamming"`, `"ppc"`
#' (positional proximity coefficient), and `"aprd"` (absolute pairwise
#' rank differences).
#' @param reverse a logical indicating if the orders should also be checked in
#' reverse order and the best value (highest correlation, lowest distance) is
#' reported. This only affect ranking-based measures and not precedence
#' invariant measures (e.g., ppc, aprd).
#' @param test a logical indicating if a correlation test should be performed.
#' @param ... Further arguments passed on to the method.
#' @return
#' - `ser_dist()` returns an object of class [dist].
#' - `ser_align()` returns a new list with elements of class
#'   [ser_permutation].
#' @author Michael Hahsler
#' @references P. Diaconis (1988): Group Representations in Probability and
#' Statistics. Institute of Mathematical Statistics, Hayward, CA.
#'
#' J.Y. Goulermas, A. Kostopoulos, and T. Mu (2016): A New Measure for
#' Analyzing and Fusing Sequences of Objects. _IEEE Transactions on
#' Pattern Analysis and Machine Intelligence_ **38**(5):833-48.
#' \doi{10.1109/TPAMI.2015.2470671}
#' @keywords cluster
#' @examples
#' set.seed(1234)
#' ## seriate dist of 50 flowers from the iris data set
#' data("iris")
#' x <- as.matrix(iris[-5])
#' x <- x[sample(1:nrow(x), 50), ]
#' rownames(x) <- 1:50
#' d <- dist(x)
#'
#' ## Create a list of different seriations
#' methods <- c("HC_single", "HC_complete", "OLO", "GW", "R2E", "VAT",
#'   "TSP", "Spectral", "SPIN", "MDS", "Identity", "Random")
#'
#' os <- sapply(methods, function(m) {
#'   cat("Doing", m, "... ")
#'   tm <- system.time(o <- seriate(d, method = m))
#'   cat("took", tm[3],"s.\n")
#'   o
#' })
#'
#' ## Compare the methods using distances. Default is based on
#' ## Spearman's rank correlation coefficient where reverse orders are
#' ## also considered.
#' ds <- ser_dist(os)
#' hmap(ds, margin = c(7,7))
#'
#' ## Compare using correlation between orders. Reversed orders have
#' ## negative correlation!
#' cs <- ser_cor(os, reverse = FALSE)
#' hmap(cs, margin = c(7,7))
#'
#' ## Compare orders by allowing orders to be reversed.
#' ## Now all but random and identity are highly positive correlated
#' cs2 <- ser_cor(os, reverse = TRUE)
#' hmap(cs2, margin=c(7,7))
#'
#' ## A better approach is to align the direction of the orders first
#' ## and then calculate correlation.
#' os_aligned <- ser_align(os)
#' cs3 <- ser_cor(os_aligned, reverse = FALSE)
#' hmap(cs3, margin = c(7,7))
#'
#' ## Compare the orders using clustering. We use Spearman's foot rule
#' ## (Manhattan distance of ranks). In order to use rank-based method,
#' ## we align the direction of the orders.
#' os_aligned <- ser_align(os)
#' ds <- ser_dist(os_aligned, method = "manhattan")
#' plot(hclust(ds))
#' @export
ser_dist <-
function(x,
  y = NULL,
  method = "spearman",
  reverse = TRUE,
  ...) {
  method <- match.arg(tolower(method), .dist_methods)

  ## make sure everything is a permutation vector
  if (!is.null(y))
    x <- list(ser_permutation_vector(x), ser_permutation_vector(y))
  else
    x <- lapply(x, ser_permutation_vector)

  if (!reverse)
    switch(
      method,
      spearman = stats::as.dist(1 - ser_cor(
        x, method = "spearman", reverse = FALSE
      )),
      kendall = stats::as.dist(1 - ser_cor(
        x, method = "kendal", reverse = FALSE
      )),

      ### Manhattan == Spearman's footrule
      manhattan = stats::dist(t(.lget_rank(x)), method = "manhattan"),
      euclidean = stats::dist(t(.lget_rank(x)), method = "euclidean"),
      hamming   = .dist_hamming(t(.lget_rank(x))),
      ppc       = as.dist(1 - ser_cor(
        x, method = "ppc", reverse = FALSE
      )),
      aprd      = stats::as.dist(.aprd(x, ...))
    )

  else
    switch(
      method,
      spearman  = stats::as.dist(1 - ser_cor(
        x, method = "spearman", reverse = TRUE
      )),
      kendall   =  stats::as.dist(1 - ser_cor(
        x, method = "kendal", reverse = TRUE
      )),

      ### Manhattan == Spearman's footrule
      manhattan = .find_best(dist(t(
        .lget_rank(.add_rev(x))
      ),
        method = "manhattan")),
      euclidean = .find_best(dist(t(
        .lget_rank(.add_rev(x))
      ),
        method = "euclidean")),
      hamming   = .find_best(.dist_hamming(t(
        .lget_rank(.add_rev(x))
      ))),

      ### positional proximity coefficient is direction invariant
      ppc       = stats::as.dist(1 - ser_cor(
        x, method = "ppc", reverse = FALSE
      )),
      aprd      = stats::as.dist(.aprd(x, ...))
    )
}

#' @rdname ser_dist
#' @export
ser_cor <- function(x,
  y = NULL,
  method = "spearman",
  reverse = TRUE,
  test = FALSE) {
  ## Note: not all .dist_methods are implemented!
  method <- match.arg(tolower(method), .dist_methods)

  ## make sure everything is a permutation vector
  if (!is.null(y))
    x <- list(ser_permutation_vector(x), ser_permutation_vector(y))
  else
    x <- lapply(x, ser_permutation_vector)

  m <- .lget_rank(x)

  if (method == "ppc") {
    if (test)
      stop("No test for association available for PPC!")
    return(.ppc(x))
  }

  ## cor based methods
  co <- stats::cor(m, method = method)
  if (reverse)
    co <- abs(co)

  ## add a correlation test?
  if (test) {
    p <- outer(1:ncol(m), 1:ncol(m), FUN =
        Vectorize(function(i, j)
          stats::cor.test(m[, i], m[, j], method = method)$p.value))
    dimnames(p) <- dimnames(co)
    attr(co, "p-value") <- p
  }

  co
}

#' @rdname ser_dist
#' @export
ser_align <- function(x, method = "spearman") {
  if (!is.list(x))
    stop("x needs to be a list with elements of type 'ser_permutation_vector'")

  x <- lapply(x, ser_permutation_vector)

  .do_rev(x, .alignment(x, method = method))
}

.dist_hamming <- function(x) {
  n <- nrow(x)
  m <- matrix(nrow = n, ncol = n)
  for (i in seq_len(n))
    for (j in seq(i, n))
      m[j, i] <- m[i, j] <- sum(x[i, ] != x[j, ])
  mode(m) <- "numeric"
  dimnames(m) <- list(rownames(x), rownames(x))
  stats::as.dist(m)
}

### make a permutation list into a rank matrix (cols are permutations)
.lget_rank <- function(x)
  sapply(x, get_rank)

### add reversed permutations to a list of permutations
.add_rev <- function(x) {
  os <- append(x, lapply(x, rev))
  names(os) <- c(labels(x), paste(labels(x), "_rev", sep = ""))
  os
}

### reverses permutations in the list given a logical indicator vector
.do_rev <- function(x, rev) {
  for (i in which(rev))
    x[[i]] <- rev(x[[i]])
  x
}

### finds the smallest distance in lists with reversed orders present
.find_best <- function(d) {
  ### find smallest values
  m <- as.matrix(d)
  n <- nrow(m) / 2
  m1 <- m[1:n, 1:n]
  m2 <- m[(n + 1):(2 * n), (n + 1):(2 * n)]
  m3 <- m[1:n, (n + 1):(2 * n)]
  m4 <- m[(n + 1):(2 * n), 1:n]
  stats::as.dist(pmin(m1, m2, m3, m4))
}

### find largest values in matrix
.find_best_max <- function(d) {
  m <- as.matrix(d)
  n <- nrow(m) / 2
  m1 <- m[1:n, 1:n]
  m2 <- m[(n + 1):(2 * n), (n + 1):(2 * n)]
  m3 <- m[1:n, (n + 1):(2 * n)]
  m4 <- m[(n + 1):(2 * n), 1:n]
  pmax(m1, m2, m3, m4)
}

### x needs to be a list of ser_permutation_vectors
### returns TRUE for sequences which should be reversed
.alignment <- function(x, method = "spearman") {
  method <- match.arg(tolower(method), .dist_methods)

  n <- length(x)

  ## calculate dist (orders + reversed orders)
  d <-
    as.matrix(ser_dist(.add_rev(x), method = method, reverse = FALSE))
  diag(d) <- NA
  for (i in 1:n) {
    d[i, n + i] <- NA
    d[n + i, i] <- NA
  }

  ## start with closest pair
  take <- which(d == min(d, na.rm = TRUE), arr.ind = TRUE)[1, ]
  #d[, c(take, (take+n) %% (2*n))] <- NA

  ## mark order and complement as taken
  d[, c(take, (take + n) %% (2 * n))] <- Inf

  ## keep adding the closest
  while (length(take) < n) {
    t2 <-
      which(d[take, ] == min(d[take, ], na.rm = TRUE), arr.ind = TRUE)[1, 2]
    #d[, c(t2,  (t2+n) %% (2*n))] <- NA

    ### closest to all
    #t2 <- which.min(colSums(d[take,], na.rm = T))
    d[, c(t2,  (t2 + n) %% (2 * n))] <- Inf

    take <- append(take, t2)
  }

  ## create indicator vector for the orders which need to be reversed
  take_ind <- logical(n)
  take_ind[take[take > n] - n] <- TRUE
  names(take_ind) <- names(x)
  take_ind
}

## Propositional Proximity Coefficient (1 - generalized corr. coef.)
## Goulermas, Kostopoulos and Mu (2016). A new measure for analyzing and fusing
## sequences of objects, IEEE Transactions on Pattern Analysis and Machine
## Intelligence 38(5):833-48.
##
## x,y ... permutation vectors (ranks)
.vppc <- Vectorize(function(x, y) {
  x <- get_rank(x)
  y <- get_rank(y)
  n <- length(x)

  #sum <- 0
  #for(j in 2:n) for(i in 1:(j-1)) sum <- sum + (x[i]-x[j])^2 * (y[i]-y[j])^2

  ## use fast matrix algebra instead
  Ax <-
    (x %*% rbind(rep_len(1, n)) - tcrossprod(cbind(rep_len(1, n)), x)) ^ 2
  Ay <-
    (y %*% rbind(rep_len(1, n)) - tcrossprod(cbind(rep_len(1, n)), y)) ^ 2
  ## note: Ay is symetric
  sum <- sum(diag(Ax %*% Ay))

  ## scale by theoretical maximum
  zapsmall(sum / (n ^ 6 / 15 - n ^ 4 / 6 + n ^ 2 / 10))
})

.ppc <- function(x)
  outer(x, x, .vppc)

# Sum of differences of rank differences
#
# distance(R, S) =
#  \sum_{i,j} | |\pi_R(i)-\pi_R(j)| - |\pi_S(i)-\pi_S(j)| |^p
#

.vaprd <- Vectorize(function(x, y, p = 2) {
  x <- get_rank(x)
  y <- get_rank(y)
  n <- length(x)

  sum <- 0
  for (j in 2:n)
    for (i in 1:(j - 1))
      sum <- sum + abs(abs(x[i] - x[j]) - abs(y[i] - y[j])) ^ p

  ## FIXME: scale by theoretical maximum?
  sum
})

.aprd <- function(x, p = 2)
  outer(x, x, .vaprd, p = p)
