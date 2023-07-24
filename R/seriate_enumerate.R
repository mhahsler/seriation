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


# utilities from package smacof
next.perm <-
  function(x)
    .C("permNext", as.double(x), as.integer(length(x)), PACKAGE = "seriation")[[1]]

are.monotone <-
  function(x, y)
    as.logical(.C(
      "isMon",
      as.double(x),
      as.double(y),
      as.integer(length(x)),
      as.integer(1),
      PACKAGE = "seriation"
    )[[4]])

.control_enumerate <- list(criterion = "Gradient_weighted",
                           verbose = FALSE)

attr(.control_enumerate, "help") <-
  list(criterion = "Criterion measure to optimize")

seriate_dist_enumerate <- function(x, control = NULL) {
  control <- .get_parameters(control, .control_enumerate)

  n <- attr(x, "Size")
  perm <- seq(n)

  best_perm <- perm
  best_crit <- Inf

  suppressWarnings(m <- as.integer(factorial(n)))
  if (is.na(m))
    stop("Number of permutations is too large.")
  k <- 0L

  if (control$verbose)
    cat("Permutation - of", m)

  repeat {
    k <- k + 1L
    if (control$verbose) {
      cat("\rPermutation", k, "of", m)
    }

    crit <-
      criterion(x,
                perm,
                method = control$criterion,
                force_loss = TRUE)

    if (crit < best_crit) {
      best_crit <- crit
      best_perm <- perm
    }

    #if (prod(perm==(n:1))==1) break
    if (k >= m)
      break
    perm <- next.perm(perm)
  }

  if (control$verbose)
    cat("\n")

  names(best_perm) <- attr(x, "Labels")[best_perm]
  best_perm
}


set_seriation_method(
  "dist",
  "Enumerate",
  seriate_dist_enumerate,
  "Enumerate all permutations",
  control = .control_enumerate,
  optimizes = .opt (NA, "set via control criterion)")
)
