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

#' Register a Genetic Algorithm Seriation Method
#'
#' Register a GA-based seriation metaheuristic for use with [seriate()].
#'
#' Registers the method `"GA"` for [seriate()]. This method can be used
#' to optimize any criterion in package \pkg{seriation}.
#'
#' The GA uses by default the ordered cross-over (OX) operator. For mutation,
#' the GA uses a mixture of simple insertion and simple inversion operators.
#' This mixed operator is created using
#' `seriation::gaperm_mixedMutation(ismProb = .8)`, where `ismProb`
#' is the probability that the simple insertion mutation operator is used. See
#' package \pkg{GA} for a description of other available cross-over and
#' mutation operators for permutations. The appropriate operator functions in
#' \pkg{GA} start with `gaperm_`.
#'
#' We warm start the GA using `"suggestions"` given by several heuristics.
#' Set `"suggestions"` to `NA` to start with a purely random initial
#' population.
#'
#' See Example section for available control parameters.
#'
#' **Note:** Package \pkg{GA} needs to be installed.
#'
#' @aliases register_GA GA ga gaperm_mixedMutation
#' @family seriation
#' @returns Nothing.
#'
#' @author Michael Hahsler
#' @references Luca Scrucca (2013): GA: A Package for Genetic Algorithms in R.
#' _Journal of Statistical Software,_ **53**(4), 1--37. URL
#' \doi{10.18637/jss.v053.i04}.
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_GA()
#' get_seriation_method("dist", "GA")
#'
#' data(SupremeCourt)
#' d <- as.dist(SupremeCourt)
#'
#' ## optimize for linear seriation criterion (LS)
#' o <- seriate(d, "GA", criterion = "LS", verbose = TRUE)
#' pimage(d, o)
#'
#' ## Note that by default the algorithm is already seeded with a LS heuristic.
#' ## This run is no warm start (no suggestions) and increase run to 100
#' o <- seriate(d, "GA", criterion = "LS", suggestions = NA, run = 100,
#'   verbose = TRUE)
#' pimage(d, o)
#'
#' o <- seriate(d, "GA", criterion = "LS", suggestions = NA, run = 100,
#'   verbose = TRUE,  )
#'
#' pimage(d, o)
#' }
#' @export
register_GA <- function() {
  check_installed("GA")

  .ga_contr <- structure(list(
    criterion = "BAR",
    suggestions = c("TSP", "QAP_LS", "Spectral"),
    selection = GA::gaperm_lrSelection,
    crossover = GA::gaperm_oxCrossover,
    mutation = gaperm_mixedMutation(.8),
    pcrossover = .8,
    pmutation = .1,
    popSize = 100,
    maxiter = 1000,
    run = 50,
    parallel = FALSE,
    verbose = FALSE
  ), help = list(
    criterion = "criterion to be optimized",
    suggestions = "seed the population with these seriation methods",
    selection = "selection operator function",
    crossover = "crossover operator function",
    mutation = "mutation operator function",
    pcrossover = "probability for crossover",
    pmutation = "ptobability of mutations",
    popSize = "population size",
    maxiter = "maximum number of generations",
    run = "stop after run generations without improvement",
    parallel = "use multiple cores?"
  ))

  GA_helper <- function(x, control) {
    n <- attr(x, "Size")

    control <- .get_parameters(control, .ga_contr)

    if (control$verbose)
      cat("\nPreparing suggestions:",
          paste0(control$suggestions, collapse = ", "), "\n")

    if (is.na(control$suggestions[1]))
      suggestions <- NULL
    else
      suggestions <- t(sapply(control$suggestions,
        function(method)
          get_order(seriate(x, method = method))))

    if (control$verbose)
      cat("\nStarting GA\n")

    # fitness function
    f <-
      function(o)
        - criterion(x, as.integer(o), method = control$criterion, force_loss = TRUE)

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
      monitor = control$verbose,
      parallel = control$parallel,
      maxiter = control$maxiter,
      run = control$run,
      maxFitness = Inf,
      popSize = control$popSize
    )

    if (control$verbose)
      if (result@iter < control$maxiter)
        cat("\nStopped early after", control$run, "iterations with no improvement! (control option 'run')\n")

    # solution may have multiple rows! Take the first solution.
    as.integer(result@solution[1, , drop = TRUE])
  }

  set_seriation_method(
    "dist",
    "GA",
    GA_helper,
    "Use a genetic algorithm to optimize for various criteria.",
    .ga_contr,
    randomized = TRUE,
    optimizes = .opt(NA, "specified as parameter criterion"),
    verbose = TRUE
  )
}


# Generates a mutation function which mixes simMutation (simple insertion)
# with ismMutation (inversion) given the probability.

#' @rdname register_GA
#' @param ismProb probability to use [GA::gaperm_ismMutation()] (inversion) versus [GA::gaperm_simMutation()] (simple insertion).
#' @export
gaperm_mixedMutation <- function(ismProb = .8) {
  function(object, parent, ...) {
    if (runif(1) > ismProb)
      GA::gaperm_simMutation(object, parent, ...)
    else
      GA::gaperm_ismMutation(object, parent, ...)
  }
}
