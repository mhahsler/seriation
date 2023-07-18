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
  cool = 0.5,
  t_min = 1e-7,
  localsearch = "LS_insert",
  try_multiplier = 5,
  t0 = NA,
  p_initial_accept = .01,
  warmstart = "Random",
  ## use "Random" for random init.
  ## try try_multiplier x n local search steps
  verbose = FALSE
)

attr(.sa_contr, "help") <- list(
  criterion = "Criterion measure to optimize",
  cool = "cooling factor (smaller means faster cooling)",
  t_min = "stopping temperature",
  localsearch = "used local search move function",
  try_multiplier = "number of local move tries per object",
  t0 = "initial temperature (if NA then it is estimated)",
  p_initial_accept = "Probability to accept a bad move at time 0 (used for t0 estimation)",
  warmstart = "permutation or seriation method for warmstart"
  )


seriate_sa <- function(x, control = NULL) {
  param <- .get_parameters(control, .sa_contr)
  n <- attr(x, "Size")

  localsearch <- get(param$localsearch)
  if (!is.function(localsearch))
    localsearch <- get(localsearch)

  crit <- param$crit

  if (is.ser_permutation(param$warmstart)) {
    .check_dist_perm(x, order = param$warmstart)
    o <- get_order(param$warmstart)
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

  iloop <- param$try_multiplier * n

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
        localsearch(o_rand),
        method = param$criterion,
        force_loss = TRUE
      )
    })

    deltas <- (z_rand - z_new)
    deltas[deltas > 0] <- NA
    avg_delta <- stats::median(deltas, na.rm = TRUE)
    t0 <- avg_delta / log(param$p_initial_accept)
  }

  nloop <- as.integer((log(param$t_min) - log(t0)) / log(param$cool))

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
      onew <- localsearch(o)
      znew <-
        criterion(x,
                  onew,
                  method = crit,
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
  "Minimize a specified seriation measure (criterion) using simulated annealing.",
  .sa_contr,
  randomized = TRUE
)
