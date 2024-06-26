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


#' @rdname criterion
#' @export
criterion.dist <- function(x,
  order = NULL,
  method = NULL,
  force_loss = FALSE,
  ...) {
  ## check dist (most C code only works with lower-triangle version)
  if (attr(x, "Diag") || attr(x, "Upper"))
    x <- as.dist(x, diag = FALSE, upper = FALSE)
  if (!is.double(x))
    mode(x) <- "double"

  ## check order
  if (!is.null(order)) {
    if (!inherits(order, "ser_permutation"))
      order <- ser_permutation(order)
    .check_dist_perm(x, order)
  } else
    order <- ser_permutation(seq(attr(x, "Size")))



  ## get methods
  if (is.null(method))
    method <- list_criterion_methods("dist")
  method <-
    lapply(method, function(m)
      get_criterion_method("dist", m))

  crit <- sapply(method,
    function(m)
      structure(m$fun(x, order, ...), names = m$name))

  if (force_loss)
    crit <- crit * sapply(
      method,
      FUN =
        function(m)
          ifelse(m$merit, -1, 1)
    )

  crit
}

#' @export
criterion.default <- criterion.dist

## Wrapper to computing the length of the order under a distance matrix,
## e.g. a tour where the leg between the first and last city is omitted.
## that this is a (Hamilton) path.
##
## Note that this corresponds to the sum of distances along the first
## off diagonal of the ordered distance matrix.

criterion_path_length <- function(x, order = NULL, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("order_length", x, order, PACKAGE = "seriation")
}

criterion_lazy_path_length <- function(x, order = NULL, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("lazy_path_length", x, order, PACKAGE = "seriation")
}

## Least squares criterion. measures the difference between the
## dissimilarities between two elements and the rank distance
## (PermutMatrix).

criterion_least_squares <- function(x, order = NULL, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("least_squares_criterion", x, order, PACKAGE = "seriation")
}

## inertia around the diagonal (see PermutMatrix)
criterion_inertia <- function(x, order = NULL, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("inertia_criterion", x, order, PACKAGE = "seriation")
}

## anti-Robinson loss functions (Streng and Schoenfelder 1978, Chen
## 2002)
## method: 1...i, 2...s, 3...w
.ar <- function(x, order = NULL, method = 1L) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("ar", x, order, as.integer(method), PACKAGE = "seriation")
}

criterion_ar_events <- function(x, order, ...)
  .ar(x, order, 1L)

criterion_ar_deviations <- function(x, order, ...)
  .ar(x, order, 2L)

#criterion_ar_weighted <- function(x, order, ...) .ar(x, order, 3L)


.rgar_contr <- structure(
  list(
    w = NULL,
    pct = 100,
    relative = TRUE
  ),
  help = list(
    w = "window size. Default is to use a pct of 100% of n",
    pct = "specify w as a percentage of n in (0,100]",
    relative = "set to FALSE to get the GAR, i.e., the absolute number of AR events in the window."
  )
)

## w \in [2,n-1]
## or pct \in [0 and 100%]; 0 -> 2 and 100 -> n-1
criterion_rgar <-
  function(x,
    order,
    w = NULL,
    pct = 100,
    relative = TRUE,
    ...) {

    if (is.null(order))
      order <- 1:attr(x, "Size")
    else
      order <- get_order(order)

    if (is.null(w)) {
      w <- floor((length(order) - 3L) * pct / 100) + 2L
      if (w < 1)
        w <- 1
    }

    if (w < 2 ||
        w >= length(order))
      stop("Window w needs to be 2 <= w < length(order) or pct needs to be 0 < pct <= 100!")
    .Call("rgar",
      x,
      order,
      as.integer(w),
      as.integer(relative),
      PACKAGE = "seriation")
  }

.bar_contr <- structure(
  list(
    b = NULL
  ),
  help = list(
    b = "band size defaults to a band of 20% of n"
  )
)

criterion_bar <- function(x, order, b = NULL, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)

  ### we default to 1/5
  if (is.null(b))
    b <- max(1, floor(length(order) / 5))

  if (b < 1 || b >= length(order))
    stop("Band size needs to be 1 <= b < length(order)!")
  .Call("bar", x, order, as.integer(b), PACKAGE = "seriation")
}

criterion_gradient_raw <- function(x, order, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("gradient", x, order, 1L, PACKAGE = "seriation")
}

criterion_gradient_weighted <- function(x, order, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)
  .Call("gradient", x, order, 2L, PACKAGE = "seriation")
}

.A_2SUM <- function(n)
  outer(
    1:n,
    1:n,
    FUN = function(i, j)
      (i - j) ^ 2
  )

criterion_2SUM <- function(x, order, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)

  # this is sum(diag(A%*%B[o,o]))
  qap::qap.obj(.A_2SUM(attr(x, "Size")), 1 / (1 + as.matrix(x)), order)
}

### Note: We use n-abs(1-j) since QAP needs positive entries in A!
.A_LS <- function(n)
  outer(
    1:n,
    1:n,
    FUN = function(i, j)
      n - abs(i - j)
  )

criterion_LS <- function(x, order, ...) {
  if (is.null(order))
    order <- 1:attr(x, "Size")
  else
    order <- get_order(order)

  # this is sum(diag(A%*%B[o,o]))
  qap::qap.obj(.A_LS(attr(x, "Size")), as.matrix(x), order)
}

# Spearman rank correlation between distances and rank differences of the order
criterion_R_dist  <- function(x, order, ...)
  abs(stats::cor(x, stats::dist(get_rank(order), "manhattan"), method = "spearman"))

### these measures are calculated on similarity matrices
criterion_ME_dist <- function(x, order, ...)
  criterion(1 / (1 + as.matrix(x)), c(order, order), "ME")
criterion_Moore_stress_dist  <- function(x, order, ...)
  criterion(1 / (1 + as.matrix(x)), c(order, order),
    "Moore_stress")
criterion_Neumann_stress_dist  <- function(x, order, ...)
  criterion(1 / (1 + as.matrix(x)), c(order, order),
    "Neumann_stress")



### register methods
set_criterion_method("dist",
  "AR_events" ,
  criterion_ar_events,
  "Anti-Robinson events: The number of violations of the anti-Robinson form (Chen, 2002).",
  FALSE)

set_criterion_method("dist",
  "AR_deviations",
  criterion_ar_deviations,
  "Anti-Robinson deviations: The number of violations of the anti-Robinson form weighted by the deviation (Chen, 2002).",
  FALSE)
## set_criterion_method("dist", "AR_weighted", criterion_ar_weighted)

set_criterion_method("dist",
  "RGAR",
  criterion_rgar,
  "Relative generalized anti-Robinson events: Counts Anti-Robinson events in a variable band of size w around the main diagonal and normalizes by the maximum of possible events (Tien et al, 2008).",
  FALSE, control = .rgar_contr)

set_criterion_method("dist", "BAR", criterion_bar,
  "Banded Anti-Robinson form criterion: Measure for closeness to the anti-Robinson form in a band of size b (Earle and Hurley, 2015).",
  FALSE,
  control = .bar_contr)

set_criterion_method("dist",
  "Gradient_raw" ,
  criterion_gradient_raw,
  "Gradient measure: Evaluates how well distances increase when moving away from the diagonal of the distance matrix (Hubert et al, 2001).",
  TRUE)

set_criterion_method(
  "dist",
  "Gradient_weighted",
  criterion_gradient_weighted,
  "Gradient measure (weighted): Evaluates how well distances increase when moving away from the diagonal of the distance matrix (Hubert et al, 2001).",
  TRUE
)

set_criterion_method("dist",
  "Path_length",
  criterion_path_length,
  "Hamiltonian path length: Sum of distances by following the permutation (Caraux and Pinloche, 2005).",
  FALSE)

set_criterion_method("dist",
  "Lazy_path_length",
  criterion_lazy_path_length,
  "Lazy path length: A weighted version of the Hamiltonian path criterion where later distances are less important (Earl and Hurley, 2015).",
  FALSE)

set_criterion_method("dist", "Inertia", criterion_inertia,
  "Inertia criterion: Measures the moment of the inertia of dissimilarity values around the diagonal of the distance matrix (Caraux and Pinloche, 2005).",
  TRUE)

set_criterion_method("dist",
  "Least_squares",
  criterion_least_squares,
  "Least squares criterion: The sum of squared differences between distances and the rank differences (Caraux and Pinloche, 2005).",
  FALSE)

set_criterion_method("dist",
  "ME",
  criterion_ME_dist,
  "Measure of effectiveness applied to the reordered similarity matrix (McCormick, 1972).",
  TRUE)

set_criterion_method("dist",
  "Rho",
  criterion_R_dist,
  "Absolute value of the Spearman rank correlation between original distances and rank differences of the order.",
  TRUE)

set_criterion_method(
  "dist",
  "Moore_stress",
  criterion_Moore_stress_dist,
  "Stress criterion (Moore neighborhood) applied to the reordered similarity matrix (Niermann, 2005).",
  FALSE
)

set_criterion_method(
  "dist",
  "Neumann_stress",
  criterion_Neumann_stress_dist,
  "Stress criterion (Neumann neighborhood) applied to the reordered similarity matrix (Niermann, 2005).",
  FALSE
)

set_criterion_method("dist",
  "2SUM",
  criterion_2SUM,
  "2-Sum Criterion: The 2-Sum loss criterion multiplies the similarity between objects with the squared rank differences (Barnard, Pothen and Simon, 1993).",
  FALSE)

set_criterion_method("dist",
  "LS",
  criterion_LS,
  "Linear Seriation Criterion: Weights the distances with the absolute rank differences (Hubert and Schultz, 1976).",
  FALSE)
