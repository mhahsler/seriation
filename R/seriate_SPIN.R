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


## SPIN (Tsafrir et al. 2005)

## Weight matrix
## pimage(create_x(n=150, sigma=20, verbose=TRUE))
create_W <- function(n, sigma, verbose = FALSE) {
  w <- function(i, j, n, sigma)
    exp(-1 * (i - j) ^ 2 / n / sigma)
  W <- outer(1:n,
    1:n,
    FUN = w,
    n = n,
    sigma = sigma)

  ## make doubly stochastic
  for (i in 1:1000) {
    #cat(i, ".")
    W <- sweep(W, MARGIN = 1, STATS = rowSums(W), "/")
    W <- sweep(W, MARGIN = 2, STATS = colSums(W), "/")
    if (all(round(rowSums(W), 5) == 1) &&
        all(round(colSums(W), 5) == 1))
      break
  }

  if (verbose)
    cat("It took", i, "iterations to make W doubly stochastic!\n")
  if (i > 999)
    warning("Weight matrix did not converge to doubly stochastic in 1000 itermation!")
  W
}

.spin_contr <- structure(
  list(
    sigma = floor(seq(20, 1, length.out = 10)),
    step = 5,
    W_function = NULL,
    verbose = FALSE
  ),
  help = list(
    sigma =  "emphasize local (small alpha) or global (large alpha) structure.",
    step = "iterations to run for each sigma value.",
    W_function = "custom function to create the weight matrix W"
  )
)

## SPIN: Neighborhood algorithms
seriate_dist_SPIN <- function(x, control = NULL) {
  param <- .get_parameters(control, .spin_contr)

  W_function <-
    if (is.null(param$W_function))
      create_W
  else
    param$W_function
  sigma <- param$sigma
  step <- param$step
  verbose <- param$verbose

  D <- as.matrix(x)
  n <- nrow(D)

  ## weight matrix
  W <- W_orig <- W_function(n, sigma[1], verbose)

  energy_best <- Inf

  for (i in 1:(length(sigma) * step)) {
    if (verbose)
      cat("Iteration", i, "... ")

    M <- D %*% W

    ## heuristic for the linear assignment problem
    ## (second argument to order breaks ties randomly)
    P <- permutation_vector2matrix(order(apply(M, MARGIN = 1, which.min), sample(1:n)))
    #if(verbose) print(table(apply(M, MARGIN = 1, which.min)))

    energy_new <- sum(diag(P %*% M))
    if (verbose)
      cat("best energy:", energy_best,
        "new energy: ", energy_new, "\n")

    ## was energy improved?
    if (energy_new < energy_best) {
      energy_best <- energy_new
      P_best <- P
    }

    ## adapt sigma
    if (!(i %% step) && i != length(sigma) * step) {
      s <- sigma[i / step + 1]
      if (verbose)
        cat("\nReducing sigma to:", s, "\n")

      W_orig <- W_function(n, s, verbose)


      ## recalculate best energy
      W <- crossprod(P, W_orig) ### t(P) %*% W
      M <- D %*% W
      energy_best <- sum(diag(P %*% M))
      if (verbose)
        cat("best energy is now:", energy_best, "\n\n")
    } else {
      W <- crossprod(P, W_orig) ### t(P) %*% W
    }
  }

  if (verbose)
    cat("Best Energy:", energy_best, "\n")
  o <- permutation_matrix2vector(P_best)
  o
}

## SPIN: Side-to-Side algorithm

## this is the weight: pimage(tcrossprod(1:n - (n+1)/2))
.spin_sts_contr <- structure(
  list(
    step = 25L,
    nstart = 10L,
    X = function(n)
      seq(n) - (n + 1) / 2,
    verbose = FALSE
  ),
  help = list(step = "iterations to run",
              nstart = "number of random restarts",
              X = "matrix to calculate the W matrix")
)

seriate_dist_SPIN_STS <- function(x, control = NULL) {
  param <- .get_parameters(control, .spin_sts_contr)

  step <- param$step
  verbose <- param$verbose
  nstart <- param$nstart
  X <- param$X

  D <- as.matrix(x)
  n <- nrow(D)

  ## X for weights W = X %*% t(X) (column vector)
  if (is.function(X))
    X <- X(n)
  if (!is.numeric(X) ||
      length(X) != n)
    stop("Invalid weight vector X.")
  W <- tcrossprod(X) ## X %*% t(X)

  .STS_run <- function() {
    if (verbose)
      cat("\nStarting new run\n")

    ## start with random permutation
    o_best <- o <- sample(1:n)
    #P_best <- P <- permutation_vector2matrix(o)
    #X_current <- crossprod(P, X)
    X_current <- X[o]
    #energy_best <- sum(diag(P %*% D %*% t(P) %*% W))
    energy_best <- sum(diag(D[o, o] %*% W))

    for (i in 1:step) {
      if (verbose)
        cat("Iteration", i, "... ")

      ## permutation matrix that orders S in descending order (break ties)
      S <- D %*% X_current
      o <- order(S, sample(1:n), decreasing = TRUE)
      #P <- permutation_vector2matrix(o)
      #X_current <- crossprod(P, X) ## t(P) %*% X
      X_current <- X[o] ## t(P) %*% X

      ## calculate energy F(P)
      #energy_new <- sum(diag(P %*% D %*% t(P) %*% W))
      energy_new <- sum(diag(D[o, o] %*% W))
      if (verbose)
        cat("best energy:", energy_best,
          "new energy: ", energy_new)

      ## was energy improved?
      if (energy_new < energy_best) {
        energy_best <- energy_new
        #P_best <- P
        o_best <- o
        if (verbose)
          cat(" - update")
      }

      if (verbose)
        cat("\n")

    }

    if (verbose)
      cat("Best Energy:", energy_best, "\n")

    #o <- permutation_matrix2vector(P_best)
    o <- o_best
    attr(o, "energy") <- energy_best
    o
  }

  res <- replicate(nstart, .STS_run(), simplify = FALSE)
  energy <- sapply(res, attr, "energy")

  if (verbose)
    cat("Overall best Energy:", min(energy), "\n")
  o <- res[[which.min(energy)]]
  o
}

set_seriation_method(
  "dist",
  "SPIN_NH",
  seriate_dist_SPIN,
  "Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005). Neighborhood algorithm to concentrate low distance values around the diagonal.",
  .spin_contr,
  optimizes = .opt(NA, "Energy")
)
set_seriation_method(
  "dist",
  "SPIN_STS",
  seriate_dist_SPIN_STS,
  "Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005). Side-to-Side algorithm which tries to push out large distance values.",
  .spin_sts_contr,
  optimizes = .opt(NA, "Energy")
)
