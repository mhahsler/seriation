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

.sgd_contr <- structure(
  list(
    criterion = "Gradient_raw",
    init = "Spectral",
    max_iter = NULL,
    localsearch = "LS_insert",
    verbose = FALSE
  ),
  help = list(
    criterion = "Criterion measure to optimize",
    init = "Start permutation or name of a seriation method",
    max_iter = "number of iterations",
    localsearch = "used local search move function"
  )
)

seriate_sgd <- function(x, control = NULL) {
  param <- .get_parameters(control, .sgd_contr)
  n <- attr(x, "Size")

  if (is.numeric(param$init)) {
    .check_dist_perm(x, order = param$init)
    o <- get_order(param$init)
  } else{
    if (param$verbose)
      cat("Obtaining initial solution via:",
          param$init, "\n")
    o <- get_order(seriate(x, method = param$init))
  }

  localsearch <- get(param$localsearch)
  if (!is.function(localsearch))
    localsearch <- get(localsearch)

  crit <- param$criterion

  max_iter <- control$max_iter
  if (is.null(max_iter))
    max_iter <- 100 * n

  z <- criterion(x, o, method = crit, force_loss = TRUE)

  if (param$verbose) {
    cat("Initial z =", z,
        "(minimize)\n")

    cat("\nTry\n")
  }

  zbest <- z

  for (i in seq(max_iter)) {
    o_new <- localsearch(o)
    z_new <-
      criterion(x,
                o_new,
                method = crit,
                force_loss = TRUE)
    delta <- z - z_new

    # we minimize, delta < 0 is a bad move
    if (delta > 0) {
      o <- o_new
      z <- z_new
      if (param$verbose)
        cat(i, "/", max_iter, "\tz =", z, "\n")
    }
  }

  o
}

set_seriation_method(
  "dist",
  "SGD",
  seriate_sgd,
  "Improve an existing solution using stochastic gradient descent.",
  .sgd_contr,
  randomized = TRUE
)
