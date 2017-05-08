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

## Brusco: simulated annealing for the Linear Seriation Criterion
seriate_dist_arsa <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    cool = 0.5,                  ## Brusco: 0.95
    tmin = 0.0001,               ## Brusco: 0.0001
    swap_to_inversion = .5,      ## Brusco: .5
    try_multiplier = 100,        ## Brusco: 100
    reps = 1L,                   ## Brusco: 20
    verbose = FALSE
  ))

  A <- as.matrix(x)
  # SUBROUTINE arsa(N, A, COOL, TMIN, NREPS, IPERM, R1, R2, D, U,
  #      S, T, SB, ZBEST, verbose)
  N <- ncol(A)
  NREPS <- as.integer(param$reps)
  IPERM <- integer(N)
#  R1 <- double(N*N/2)
#  R2 <- double(N*N/2)
  D <- double(N*N)
  U <- integer(N)
  S <- integer(N)
  T <- integer(NREPS*N)
  SB <- integer(N)
  ZBEST <- double(1)

  ret <- .Fortran("arsa", N, A,
    as.numeric(param$cool), as.numeric(param$tmin), NREPS,
    IPERM, D, U, S, T, SB, ZBEST, as.numeric(param$swap_to_insertion),
    as.numeric(param$try_multiplier), as.logical(param$verbose),
    PACKAGE="seriation")

  o <- ret[[6]]
  names(o) <- labels(x)[o]

  ### ARSA returns all 0's in some cases
  if(all(o == 0)) {
    o <- 1:N
    warning("ARSA has returned an invalid permutation vector! Check the supplied dissimilarity matrix.")
  }

  o
}


## Brusco: branch-and-bound - unweighted row gradient
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


## Brusco: branch-and-bound - weighted row gradient
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

set_seriation_method("dist", "ARSA", seriate_dist_arsa,
  "Minimize the linear seriation criterion using simulated annealing")
set_seriation_method("dist", "BBURCG", seriate_dist_bburcg,
  "Minimize the unweighted row/column gradient by branch-and-bound")
set_seriation_method("dist", "BBWRCG", seriate_dist_bbwrcg,
  "Minimize the weighted row/column gradient by branch-and-bound")
