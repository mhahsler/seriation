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


## Simulated annealing reimplimentation following 'arsa.f' by Brusco et al.
## can use any criterion function

#' Neighborhood functions for Seriation Method SA
#'
#' Definition of different local neighborhood functions for the method `"SA"` for [seriate()].
#'
#' Local neighborhood functions are `LS_insert`, `LS_swap`, `LS_reverse`, and `LS_mix`
#'  (1/3 insertion, 1/3 swap and 1/3 reverse). Any neighborhood function can be defined.
#' @name LS
#' @aliases LS
#' @param o an integer vector with the order
#' @param pos random positions used for the local move.
#' @returns returns the new order vector representing the random neighbor.
NULL

#' @rdname LS
#' @export
LS_swap <- function(o, pos = sample.int(length(o), 2)) {
  tmp <- o[pos[1]]
  o[pos[1]] <- o[pos[2]]
  o[pos[2]] <- tmp
  o
}

### insert pos[1] in pos[2]
#' @rdname LS
#' @export
LS_insert <- function(o, pos = sample.int(length(o), 2)) {
  append(o[-pos[1]], o[pos[1]], after = pos[2] - 1)
}

#' @rdname LS
#' @export
LS_reverse <- function(o, pos = sample.int(length(o), 2)) {
  o[pos[1]:pos[2]] <- o[pos[2]:pos[1]]
  o
}

#' @rdname LS
#' @export
LS_mixed <- function(o, pos = sample.int(length(o), 2)) {
  switch(sample.int(3, 1),
         LS_swap(o, pos),
         LS_insert(o, pos),
         LS_reverse(o, pos))
}

.sa_contr <- list(
  criterion = "Gradient_raw",
  warmstart = "Spectral",
  ## use "Random" for random init.
  localsearch = LS_insert,
  cool = 0.5,
  tmin = 1e-7,
  nlocal = 5,
  t0 = NA,
  pinitialaccept = .01,
  ## try nlocal x n local search steps
  verbose = FALSE
)

seriate_sa <- function(x, control = NULL) {
  param <- .get_parameters(control, .sa_contr)
  n <- attr(x, "Size")

  if (is.na(param$warmstart))
    param$warmstart <- "RANDOM"

  if (is.numeric(param$init)) {
    .check_dist_perm(x, order = param$init)
  } else{
    if (param$verbose)
      cat("Obtaining initial solution via:",
          param$warmstart, "\n")
    o <- get_order(seriate(x, method = param$warmstart))
  }

  z <- criterion(x, o, method = param$criterion, force_loss = TRUE)
  if (param$verbose)
    cat("Initial z =", z,
        "(minimize)\n")

  iloop <- param$nlocal * n

  t0 <- param$t0
  if (is.na(t0)) {
    # find the starting temperature. Set the probability of the average
    # (we use median) uphill move to pinitaccept.
    o_rand <- sample(n)
    z_rand <-  criterion(x,
                         o_rand,
                         method = param$criterion,
                         force_loss = TRUE)
    z_new <- replicate(iloop, expr = {
      criterion(
        x,
        param$localsearch(o_rand),
        method = param$criterion,
        force_loss = TRUE
      )
    })

    deltas <- (z_rand - z_new)
    deltas[deltas > 0] <- NA
    avg_delta <- stats::median(deltas, na.rm = TRUE)
    t0 <- avg_delta / log(param$pinitialaccept)
  }

  nloop <- as.integer((log(param$tmin) - log(t0)) / log(param$cool))

  if (t0 <= 0) {
    t0 <- 0
    nloop <- 1L
  }

  if (param$verbose)
    cat("Use t0 =",
        t0,
        "resulting in",
        nloop,
        "iterations with",
        iloop,
        "tries each\n\n")

  zbest <- z
  temp <- t0

  for (i in 1:nloop) {
    m <- 0L

    for (j in 1:iloop) {
      onew <- param$localsearch(o)
      znew <-
        criterion(x,
                  onew,
                  method = param$criterion,
                  force_loss = TRUE)
      delta <- z - znew

      # we minimize, delta < 0 is a bad move
      if (delta > 0 || temp > 0 && runif(1) < exp(delta / temp)) {
        o <- onew
        z <- znew
        m <- m + 1L
      }
    }

    if (param$verbose) {
      cat(
        i,
        "/",
        nloop,
        "\ttemp =",
        signif(temp, 3),
        "\tz =",
        z,
        "\t accepted moves =",
        m,
        "/",
        iloop,
        "\n"
      )
    }

    temp <- temp * param$cool
  }

  o
}

set_seriation_method(
  "dist",
  "GSA",
  seriate_sa,
  paste0(
    "Minimize a specified seriation measure (criterion) using simulated annealing.\n",
    "Control parameters:\n",
    " - criterion to optimize\n",
    " - init (initial order; use \"Random\" for no warm start\n",
    " - localsearch (neighborhood function; Built-in functions are LS_insert, LS_swap, LS_reverse, and LS_mix (1/3 insertion, 1/3 swap and 1/3 reverse)\n",
    " - cool (cooling rate)\n",
    " - tmin (minimum temperature)\n",
    " - swap_to_inversion (proportion of swaps to inversions)\n",
    " - nlocal (number of objects times nlocal is the number of search tries per temperature\n"
  ),
  .sa_contr,
  randomized = TRUE
)
