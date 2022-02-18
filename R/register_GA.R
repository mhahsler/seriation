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



#' Register a Genetic Algorithm Seriation Method
#'
#' Register a GA-based seriation metaheuristic for use with [seriate()].
#'
#' Registers the method \code{"GA"} for [seriate()]. This method can be used
#' to optimize any criterion in package \pkg{seriation}.
#'
#' \code{control} for
#' \code{seriate} with method \code{"GA"} accepts the following parameters:
#' - "criterion" criterion to optimize. Default: BAR
#' - "suggestions" suggestions to warm start the GA.
#'   \code{NA} means no warm start. Default: TSP, QAP_LS and Spectral.
#' - "selection" Selection operator.
#' - "crossover" Crossover operator.
#' - "mutation" Mutation operator. Default: a
#'   mixture of the simple insertion (80% chance) and simple inversion (20%
#'   chance) operators.
#' - "pmutation" probability for permutations. Default: .5
#' - "pcrossover" probability for crossover. Default: .2
#' - "popsize" the population size. Default: 100
#' - "maxiter" maximum number of generations. Default: 1000
#' - "run" stop after \code{run} generations without improvement. Default: 50
#' - "parallel" use multiple cores? Default: TRUE
#' - "verbose" a logical; report progress? Default: TRUE
#'
#' The GA uses by default the ordered cross-over (OX) operator. For mutation,
#' the GA uses a mixture of simple insertion and simple inversion operators.
#' This mixed operator is created using
#' \code{seriation::gaperm_mixedMutation(ismProb = .8)}, where \code{ismProb}
#' is the probability that the simple insertion mutation operator is used. See
#' package \pkg{GA} for a description of other available cross-over and
#' mutation operators for permutations. The appropriate operator functions in
#' \pkg{GA} start with \code{gaperm_}.
#'
#' We warm start the GA using \code{"suggestions"} given by several heuristics.
#' Set \code{"suggestions"} to \code{NA} to start with a purely random initial
#' population.
#'
#' \bold{Note:} Package \pkg{GA} needs to be installed.
#'
#' @aliases register_GA GA ga gaperm_mixedMutation
#' @family seriation
#' @returns Nothing.
#'
#' @author Michael Hahsler
#' @references Luca Scrucca (2013): GA: A Package for Genetic Algorithms in R.
#' \emph{Journal of Statistical Software,} \bold{53}(4), 1--37. URL
#' \doi{10.18637/jss.v053.i04}.
#' @keywords optimize cluster
#' @examples
#'
#' \dontrun{
#' register_GA()
#' get_seriation_method("dist", "GA")
#'
#' d <- dist(random.robinson(50, pre=TRUE, noise=.1))
#'
#' ## use default settings: Banded AR form
#' o <- seriate(d, "GA")
#' pimage(d, o)
#'
#' ## optimize for linear sertiation criterion (LS)
#' o <- seriate(d, "GA", control = list(criterion = "LS"))
#' pimage(d, o)
#'
#' ## no warm start
#' o <- seriate(d, "GA", control = list(criterion = "LS", suggestions = NA))
#' pimage(d, o)
#' }
#' @export
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
