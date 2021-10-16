#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
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

## register GA for seriation

# Generates a mutation function which mixes simMutation (simple insertion)
# with ismMutation (inversion) given the probability.
gaperm_mixedMutation <- function(ismProb = .8) {
  function(object, parent, ...) {
    if (runif(1) > ismProb)
      GA::gaperm_simMutation(object, parent, ...)
    else
      GA::gaperm_ismMutation(object, parent, ...)
  }
}

register_GA <- function() {
  check_installed("GA")

  .ga_contr <- list(
    criterion = "BAR",
    suggestions = c("TSP", "QAP_LS", "Spectral"),
    selection = GA::gaperm_nlrSelection,
    crossover = GA::gaperm_oxCrossover,
    mutation = gaperm_mixedMutation(.8),
    pcrossover = .2,
    pmutation = .5,
    popSize = 100,
    maxiter = 1000,
    run = 50,
    parallel = TRUE,
    verbose = TRUE
  )

  GA_helper <- function(x, control) {
    n <- attr(x, "Size")

    control <- .get_parameters(control, .ga_contr)

    if (control$verbose)
      cat("\nPreparing suggestions\n")

    if (is.na(control$suggestions[1]))
      suggestions <- NULL
    else
      suggestions <- t(sapply(control$suggestions,
        function(method)
          get_order(seriate(x, method = method))))

    if (control$verbose)
      cat("\nStarting GA\n")

    ### FIXME: need to be able to set bandwidth for BAR
    # fitness function
    f <-
      function(o)
        - criterion(x, o, method = control$criterion, force_loss = TRUE)

    result <- GA::ga(
      type = "permutation",
      fitness = f,
      lower = rep(1L, times = n),
      upper = rep(n, times = n),
      selection = control$selection,
      mutation = control$mutation,
      crossover = control$crossover,
      pmutation = control$pmutation,
      pcrossover = control$pcrossover,
      suggestions = suggestions,
      names = as.character(1:n),
      monitor = if (control$verbose)
        GA::gaMonitor
      else
        NULL,
      parallel = control$parallel,
      maxiter = control$maxiter,
      run = control$run,
      maxFitness = Inf,
      popSize = control$popSize
    )

    as.integer(result@solution[1,])
  }

  set_seriation_method(
    "dist",
    "GA",
    GA_helper,
    "Use a genetic algorithm to optimize for various criteria.",
    .ga_contr
  )
}
