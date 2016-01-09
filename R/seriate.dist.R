#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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



## seriate dist objects

seriate.dist <- function(x, method = "ARSA", control = NULL, ...) {
    if(!all(x>=0)) stop("Negative distances not supported!")

    ## add ... to control
    control <- c(control, list(...))

    ## check x
    if(any(is.na(x))) stop("NAs not allowed in x!")
    if(any(x<0)) stop("No negative values allowed in x!")

    if(!is.character(method) || (length(method) != 1L))
      stop("Argument 'method' must be a character string.")
    method <- get_seriation_method("dist", method)

    if(!is.null(control$verbose) && control$verbose) cat(method$name, ": ",
      method$description, "\n", sep="")

    order <- method$fun(x, control = control)

    ser_permutation(ser_permutation_vector(order, method = method$name))
  }

## uses a sequence of correlation matrices and finds  the first matrix
## with rank 2. The elements are projected into the plane spanned by the
## first two eigenvectors. All points are lying on a ellipse. The order
## of the elements on the ellipse is returned (see Chen 2002).
seriate_dist_chen <- function(x, control = NULL){
  .get_parameters(control, NULL)

  x <- as.matrix(x)

  rank <- qr(x)$rank

  ## find the first correlation matrix of rank 2
  n <- 0
  while(rank > 2){
    x <- cor(x)
    n <- n + 1
    rank <- qr(x)$rank
  }

  ## project the matrix on the first 2 eigenvectors
  e <- eigen(x)$vectors[,1:2]

  ## extract the order
  ## chen says that he uses the one of the two possible cuts
  ## that separate the points at rank 1. Since the points just
  ## separate further towards right and left, cutting on the vertical
  ## axis of the ellipse yields the same result.

  right <- which(e[,1] >= 0)
  right <- right[order(e[right,2], decreasing = TRUE)]
  left <- which(e[,1] < 0)
  left <- left[order(e[left,2])]

  o <- c(right,left)
  names(o) <- labels(x)[o]
  o
}


## Bridge to package tsp
seriate_dist_tsp <- function(x, control = NULL){
  ## add a dummy city for cutting
  tsp <- insert_dummy(TSP(x), n = 1, label = "cut_here")

  if(is.null(control))
    control <- list(
      method="arbitrary insertion",
      rep = 10,
      two_opt = TRUE
    )

  tour <- solve_TSP(tsp, method = control$method,
    control = control$control)

  o <- cut_tour(tour, cut = "cut_here", exclude_cut = TRUE)
  names(o) <- labels(x)[o]
  o
}


## Multidimensional scaling
seriate_dist_mds <- function(x, control = NULL){
  control <- .get_parameters(control, list(
    method = "cmdscale"
  ))

  if(control$method == "cmdscale" ) {
    sc <- cmdscale(x, k=1)
    return(order(sc[,1]))

  }else if(control$method == "isoMDS"){
    sc <- MASS::isoMDS(x+1e-6, trace = FALSE, k=1)
    return(order(sc$points[,1]))

  }else if(control$method == "sammon") {
    sc <- MASS::sammon(x+1e-6, trace = FALSE, k=1)
    return(order(sc$points[,1]))

  }else stop("unknown method")

}

seriate_dist_mds_metric <- function(x, control = NULL)
  seriate_dist_mds(x, control=list(method="cmdscale"))
seriate_dist_mds_nonmetric <- function(x, control = NULL)
  seriate_dist_mds(x, control=list(method="isoMDS"))

## Angle between the first 2 PCS. Fiendly (2002)
seriate_dist_angle <- function(x, control = NULL) {
  .get_parameters(control, NULL)

  sc <- cmdscale(x, k=2)
  .order_angle(sc)
}


## Hierarchical clustering related seriations
.hclust_helper <- function(d, control = NULL){
  control <- .get_parameters(control, list(
    hclust = NULL,
    method = "average"
    ))

  if(!is.null(control$hclust)) return(control$hclust)
  return(hclust(d, method = control$method))
}

seriate_dist_hc <- function(x, control = NULL) .hclust_helper(x, control)
seriate_dist_hc_single <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="single"))
seriate_dist_hc_average <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="average"))
seriate_dist_hc_complete <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="complete"))
seriate_dist_hc_ward <- function(x, control = NULL)
  .hclust_helper(x, control=list(method="ward.D2"))

## workhorses are in seriation.hclust
seriate_dist_gw <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method="GW")
seriate_dist_gw_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method="GW")
seriate_dist_gw_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method="GW")
seriate_dist_gw_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method="GW")
seriate_dist_gw_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method="GW")


seriate_dist_olo <- function(x, control = NULL)
  reorder(seriate_dist_hc(x, control), x, method="OLO")
seriate_dist_olo_single <- function(x, control = NULL)
  reorder(seriate_dist_hc_single(x, control), x, method="OLO")
seriate_dist_olo_average <- function(x, control = NULL)
  reorder(seriate_dist_hc_average(x, control), x, method="OLO")
seriate_dist_olo_complete <- function(x, control = NULL)
  reorder(seriate_dist_hc_complete(x, control), x, method="OLO")
seriate_dist_olo_ward <- function(x, control = NULL)
  reorder(seriate_dist_hc_ward(x, control), x, method="OLO")

## brusco: simulated annealing for anti-robinson
seriate_dist_arsa <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    cool = 0.5,
    tmin = 0.1,
    nreps = 1L,
    verbose = FALSE
  ))

  A <- as.matrix(x)
  # SUBROUTINE arsa(N, A, COOL, TMIN, NREPS, IPERM, R1, R2, D, U,
  #      S, T, SB, verbose)
  N <- ncol(A)
  IPERM <- integer(N)
  R1 <- double(N*N/2)
  R2 <- double(N*N/2)
  D <- double(N*N)
  U <- integer(N)
  S <- integer(N)
  T <- integer(100*N)
  SB <- integer(N)

  ret <- .Fortran("arsa", N, A, param$cool, param$tmin, param$nreps, IPERM,
    R1, R2, D, U, S, T, SB, param$verbose, PACKAGE="seriation")

  o <- ret[[6]]
  names(o) <- labels(x)[o]

  ### ARSA returns all 0's in some cases
  if(all(o == 0)) {
    o <- 1:N
    warning("ARSA has returned an invalid permutation vector! Check the supplied dissimilarity matrix.")
  }

  o
}


## brusco: branch-and-bound - unweighted row gradient
seriate_dist_bburcg <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    eps = 1e-7,
    verbose = FALSE
  ))

  A <- as.matrix(x)
  N <- ncol(A)

  # SUBROUTINE bburcg(N, A, EPS, X, Q, D, DD, S, UNSEL, IVERB)
  X <- integer(N)
  Q <- integer(N)
  D <- integer(N*N*N)
  DD <- integer(N*N*N)
  S <- integer(N)
  UNSEL <- integer(N)

  ret <- .Fortran("bburcg", N, A, param$eps, X, Q, D, DD, S, UNSEL,
    param$verbose)

  o <- ret[[4]]
  names(o) <- labels(x)[o]
  o
}


## brusco: branch-and-bound - weighted row gradient
seriate_dist_bbwrcg <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    eps = 1e-7,
    verbose = FALSE
  ))

  A <- as.matrix(x)
  N <- ncol(A)

  # SUBROUTINE bbwrcg(N, A, EPS, X, Q, D, DD, S, UNSEL, IVERB)
  X <- integer(N)
  Q <- integer(N)
  D <- double(N*N*N)
  DD <- double(N*N*N)
  S <- integer(N)
  UNSEL <- integer(N)

  ret <- .Fortran("bbwrcg", N, A, param$eps, X, Q, D, DD, S, UNSEL,
    param$verbose)

  o <- ret[[4]]
  names(o) <- labels(x)[o]
  o
}

seriate_dist_identity <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  o <- 1:attr(x, "Size")
  names(o) <- labels(x)
  o
}

seriate_dist_random <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  o <- 1:attr(x, "Size")
  names(o) <- labels(x)
  sample(o)
}

## VAT: a tool for visual assessment of (cluster) tendency
## Bezdek, J.C., Hathaway, R.J.
## Proceedings of the 2002 International Joint Conference on
## Neural Networks, 2002. IJCNN '02. (Volume:3)
seriate_dist_VAT <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  D <- as.matrix(x)
  N <- nrow(D)
  P <- rep(NA_integer_, N)
  I <- rep(FALSE, N)
  ### J is !I

  i <- which(D == max(D, na.rm = TRUE), arr.ind = TRUE)[1,1]
  P[1] <- i
  I[i] <- TRUE

  for(r in 2:N) {
    D2 <- D[I,!I, drop=FALSE]
    j <- which(D2 == min(D2, na.rm = TRUE), arr.ind = TRUE)[1,2]
    j <- which(!I)[j]
    P[r] <- j
    I[j] <- TRUE
  }

  names(P) <- labels(x)[P]
  P
}


## Spectral Seriation
## Ding, C. and Xiaofeng He (2004): Linearized cluster assignment via
## spectral orderingProceedings of the Twenty-first.
## International Conference on Machine learning (ICML '04)

## Minimizes: sum_{i,j} (i-j)^2 * d_{pi_i,pi_j}

seriate_dist_spectral <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  ### calculate Laplacian
  W <- 1/(1+as.matrix(x))
  D <- diag(rowSums(W))
  L <- D - W

  ## Fielder vector q1 is eigenvector with the smallest eigenvalue
  ## eigen reports eigenvectors/values in decreasing order
  q <- eigen(L)
  fielder <- q$vectors[,ncol(W)-1L]
  o <- order(fielder)
  names(o) <- names(x)[o]
  o
}

seriate_dist_spectral_norm <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  ### calculate normalized Laplacian
  W <- 1/(1+as.matrix(x))
  D_sqrt<- diag(rowSums(1/W^.5))
  L <- D_sqrt %*% W %*% D_sqrt

  z <- eigen(L)$vectors
  q <- D_sqrt %*% z

  ## look for the vector with the largest eigenvalue
  largest_ev <- q[,2L]
  o <- order(largest_ev)
  names(o) <- names(x)[o]
  o
}

## SPIN (Tsafrir et al. 2005)

## Weight matrix
## pimage(create_x(n=150, sigma=20, verbose=TRUE))
create_W <- function(n, sigma, verbose=FALSE) {
  w <- function(i, j, n, sigma) exp(-1*(i-j)^2/n/sigma)
  W <- outer(1:n, 1:n, FUN = w, n=n, sigma=sigma)

  ## make doubly stochastic
  for(i in 1:1000) {
    #cat(i, ".")
    W <- sweep(W, MARGIN = 1, STATS = rowSums(W), "/")
    W <- sweep(W, MARGIN = 2, STATS = colSums(W), "/")
    if(round(rowSums(W), 5) == 1 && round(colSums(W), 5) == 1) break
  }

  if(verbose) cat("It took", i, "iterations to make W doubly stochastic!\n")
  if(i >999) warning("Weight matrix did not converge to doubly stochastic in 1000 itermation!")
  W
}

## SPIN: Neighborhood algorithms
seriate_dist_SPIN <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    sigma = seq(20,1, length.out = 10),
    step = 5,
    W_function = NULL,
    verbose = FALSE
  ))

  W_function <- if(is.null(param$W_function)) create_W else param$W_function
  sigma <- param$sigma
  step <- param$step
  verbose <- param$verbose

  D <- as.matrix(x)
  n <- nrow(D)

  ## weight matrix
  W <- W_orig <- W_function(n, sigma[1], verbose)

  energy_best <- Inf

  for(i in 1:(length(sigma)*step)) {
    if(verbose) cat("Iteration", i, "... ")

    M <- D %*% W

    ## heuristic for the linear assignment problem
    ## (second argument to order breakes ties randomly)
    P <- permutation_vector2matrix(
      order(apply(M, MARGIN = 1, which.min), sample(1:n)))
    #if(verbose) print(table(apply(M, MARGIN = 1, which.min)))

    energy_new <- sum(diag(P %*% M))
    if(verbose) cat("best energy:", energy_best,
      "new energy: ", energy_new, "\n")

    ## was energy improved?
    if(energy_new < energy_best) {
      energy_best <- energy_new
      P_best <- P
    }

    ## adapt sigma
    if(!(i %% step) && i != length(sigma)*step) {
      s <- sigma[i/step+1]
      if(verbose) cat("\nReducing sigma to:", s, "\n")

      W_orig <- W_function(n, s, verbose)


      ## recalculate best energy
      W <- crossprod(P, W_orig) ### t(P) %*% W
      M <- D %*% W
      energy_best <- sum(diag(P %*% M))
      if(verbose) cat("best energy is now:", energy_best, "\n\n")
    }else {
      W <- crossprod(P, W_orig) ### t(P) %*% W
    }
  }

  if(verbose) cat("Best Energy:", energy_best, "\n")
  o <- permutation_matrix2vector(P_best)
  names(o) <- names(x)[o]
  o
}

## SPIN: Side-to-Side algorithm

## this is the weight: pimage(tcrossprod(1:n - (n+1)/2))
seriate_dist_SPIN_STS <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    step = 25,
    nstart = 10,
    X = function(n) 1:n - (n+1)/2,
    verbose = FALSE
  ))

  step <- param$step
  verbose <- param$verbose
  nstart <- param$nstart
  X <- param$X

  D <- as.matrix(x)
  n <- nrow(D)

  ## X for weights W = X %*% t(X) (colunm vector)
  if(is.function(X)) X <- X(n)
  if(!is.numeric(X) || length(X) != n) stop("Invalid weight vector X.")
  W <- tcrossprod(X) ## X %*% t(X)

  .STS_run <- function() {
    if(verbose) cat("\nStarting new run\n")

    ## start with random permutation
    o_best <- o <- sample(1:n)
    #P_best <- P <- permutation_vector2matrix(o)
    #X_current <- crossprod(P, X)
    X_current <- X[o]
    #energy_best <- sum(diag(P %*% D %*% t(P) %*% W))
    energy_best <- sum(diag(D[o,o] %*% W))

    for(i in 1:step) {
      if(verbose) cat("Iteration", i, "... ")

      ## permutation matrix that orders S in descending order (break ties)
      S <- D %*% X_current
      o <- order(S, sample(1:n), decreasing = TRUE)
      #P <- permutation_vector2matrix(o)
      #X_current <- crossprod(P, X) ## t(P) %*% X
      X_current <- X[o] ## t(P) %*% X

      ## calculate energy F(P)
      #energy_new <- sum(diag(P %*% D %*% t(P) %*% W))
      energy_new <- sum(diag(D[o,o] %*% W))
      if(verbose) cat("best energy:", energy_best,
        "new energy: ", energy_new)

      ## was energy improved?
      if(energy_new < energy_best) {
        energy_best <- energy_new
        #P_best <- P
        o_best <- o
        if(verbose) cat(" - update")
      }

      if(verbose) cat("\n")

    }

    if(verbose) cat("Best Energy:", energy_best, "\n")

    #o <- permutation_matrix2vector(P_best)
    o <- o_best
    attr(o, "energy") <- energy_best
    o
  }

  res <- replicate(nstart, .STS_run(), simplify = FALSE)
  energy <- sapply(res, attr, "energy")

  if(verbose) cat("Overall best Energy:", min(energy), "\n")
  o <- res[[which.min(energy)]]
  names(o) <- names(x)[o]
  o
}

## QAP 2SUM seriation
seriate_dist_2SUM <- function(x, control = NULL) {
  ## param are passed on to QAP

  do.call(qap::qap, c(list(A = .A_2SUM(attr(x, "Size")),
    B = 1/(1+as.matrix(x))), control))
}

## QAP Linear seriation
seriate_dist_LS <- function(x, control = NULL) {
  ## param are passed on to QAP

  do.call(qap::qap, c(list(A = .A_LS(attr(x, "Size")),
    B = as.matrix(x)), control))
}

## QAP Inertia
seriate_dist_Inertia <- function(x, control = NULL) {
  ## param are passed on to QAP
  n <- attr(x, "Size")

  ## inertia uses the same A matrix as 2SUM
  do.call(qap::qap, c(list(A = n^2-.A_2SUM(n),
    B = as.matrix(x)), control))
}

## QAP BAR
seriate_dist_BAR <- function(x, control = NULL) {
  ## param are passed on to QAP

  .A_BAR <- function(n, b) {
    b <- floor(b)
    if(b<1 || b>=n) stop("b: needs to be 1<=b<n!")
    A <- b+1- outer(1:n,1:n, FUN = function(i,j) abs(i-j))
    A[A<0] <- 0
    A
  }

  n <- attr(x, "Size")
  if(is.null(control$b)) b <- max(1, floor(n/5))
  else {
    b <- control$b
    control$b <- NULL
  }

  ## inertia uses the same A matrix as 2SUM
  do.call(qap::qap, c(list(A = .A_BAR(n, b = b),
    B = as.matrix(x)), control))
}





set_seriation_method("dist", "Identity", seriate_dist_identity,
  "Identity permutation")
set_seriation_method("dist", "Random", seriate_dist_random,
  "Random permutation")

set_seriation_method("dist", "ARSA", seriate_dist_arsa,
  "Minimize Anti-Robinson events using simulated annealing")
set_seriation_method("dist", "BBURCG", seriate_dist_bburcg,
  "Minimize the unweighted row/column gradient by branch-and-bound")
set_seriation_method("dist", "BBWRCG", seriate_dist_bbwrcg,
  "Minimize the weighted row/column gradient by branch-and-bound")

set_seriation_method("dist", "TSP", seriate_dist_tsp,
  "Minimize Hamiltonian path length with a TSP solver")

#set_seriation_method("dist", "Chen", seriate_dist_chen,
#  "Rank-two ellipse seriation")
set_seriation_method("dist", "R2E", seriate_dist_chen,
  "Rank-two ellipse seriation")

set_seriation_method("dist", "MDS", seriate_dist_mds,
  "MDS")
set_seriation_method("dist", "MDS_metric", seriate_dist_mds_metric,
  "MDS (metric)")
set_seriation_method("dist", "MDS_nonmetric", seriate_dist_mds_nonmetric,
  "MDS (non-metric)")

set_seriation_method("dist", "MDS_angle", seriate_dist_angle,
  "MDS (angle)")

set_seriation_method("dist", "HC", seriate_dist_hc,
  "Hierarchical clustering")
set_seriation_method("dist", "HC_single", seriate_dist_hc_single,
  "Hierarchical clustering (single link)")
set_seriation_method("dist", "HC_complete", seriate_dist_hc_complete,
  "Hierarchical clustering (complete link)")
set_seriation_method("dist", "HC_average", seriate_dist_hc_average,
  "Hierarchical clustering (avg. link)")
set_seriation_method("dist", "HC_ward", seriate_dist_hc_ward,
  "Hierarchical clustering (Ward's method)")

set_seriation_method("dist", "GW", seriate_dist_gw,
  "Hierarchical clustering reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_single", seriate_dist_gw_single,
  "Hierarchical clustering (single link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_average", seriate_dist_gw_average,
  "Hierarchical clustering (avg. link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_complete", seriate_dist_gw_complete,
  "Hierarchical clustering (complete link) reordered by Gruvaeus and Wainer heuristic")
set_seriation_method("dist", "GW_ward", seriate_dist_gw_ward,
  "Hierarchical clustering (Ward's method) reordered by Gruvaeus and Wainer heuristic")


set_seriation_method("dist", "OLO", seriate_dist_olo,
  "Hierarchical clustering (single link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_single", seriate_dist_olo_single,
  "Hierarchical clustering with optimal leaf ordering")
set_seriation_method("dist", "OLO_average", seriate_dist_olo_average,
  "Hierarchical clustering (avg. link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_complete", seriate_dist_olo_complete,
  "Hierarchical clustering (complete link) with optimal leaf ordering")
set_seriation_method("dist", "OLO_ward", seriate_dist_olo_ward,
  "Hierarchical clustering (Ward's method) with optimal leaf ordering")

set_seriation_method("dist", "VAT", seriate_dist_VAT,
  "Visual assesment of clustering tendency (VAT)")

set_seriation_method("dist", "Spectral", seriate_dist_spectral,
  "Spectral seriation")
set_seriation_method("dist", "Spectral_norm", seriate_dist_spectral_norm,
  "Spectral seriation (normalized)")

set_seriation_method("dist", "SPIN_NH", seriate_dist_SPIN,
  "SPIN (Neighborhood algorithm)")
set_seriation_method("dist", "SPIN_STS", seriate_dist_SPIN_STS,
  "SPIN (Side-to-Side algorithm)")

set_seriation_method("dist", "QAP_2SUM", seriate_dist_2SUM,
  "QAP (2-SUM)")
set_seriation_method("dist", "QAP_LS", seriate_dist_LS,
  "QAP (Linear Seriation)")
set_seriation_method("dist", "QAP_BAR", seriate_dist_BAR,
  "QAP (Banded anti-Robinson form)")
set_seriation_method("dist", "QAP_Inertia", seriate_dist_Inertia,
  "QAP (Inertia)")
