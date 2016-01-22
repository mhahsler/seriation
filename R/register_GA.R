#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
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

## registers seriation methods and criteria from package DendSer

# Generates a mutation function which mixes simMutation (simple insertion)
# with ismMutation (inversion) given the probability.
gaperm_mixedMutation <- function(ismProb = .8) {
  function(object, parent, ...) {
    if(runif(1)>ismProb) GA::gaperm_simMutation(object, parent, ...)
    else GA::gaperm_ismMutation(object, parent, ...)
  }
}

register_GA <- function() {
  if(!.installed("GA")) stop("Package 'GA' needs to be  installed!")

  GA_helper <- function(x, control) {
    n <- attr(x, "Size")

    control <- .get_parameters(control, list(
      criterion = "BAR",
      suggestions = c("TSP", "QAP_LS", "Spectral"),
      selection = GA::gaperm_nlrSelection,
      crossover = GA::gaperm_oxCrossover,
      mutation = seriation::gaperm_mixedMutation(.8),
      pcrossover = .2,
      pmutation = .5,
      popSize = 100,
      maxiter = 1000,
      run = 50,
      parallel = TRUE,
      verbose = TRUE
    ))

    if(control$verbose) cat("\nPreparing suggestions\n")

    if(is.na(control$suggestions[1])) suggestions <- NULL
    else suggestions <- t(sapply(control$suggestions,
      function(method) get_order(seriate(x, method = method))))

    if(control$verbose) cat("\nStarting GA\n")

    ### FIXME: handle ...
    ## max fitness
    crit2fit <- function(x, method) {
      if(seriation::get_criterion_method("dist", method)$merit)
        function(o) criterion(x, o, method)
      else
        function(o) -criterion(x, o, method)
    }

    f <- crit2fit(x, control$criterion)

    result <- GA::ga(type="permutation",
      fitness=f,
      min=rep(1L, times = n),
      max=rep(n, times = n),
      selection = control$selection,
      mutation = control$mutation,
      crossover = control$crossover,
      pmutation = control$pmutation,
      pcrossover = control$pcrossover,
      suggestions = suggestions,
      names=as.character(1:n),
      monitor = if(control$verbose) GA::gaMonitor else NULL,
      parallel = control$parallel,
      maxiter = control$maxiter,
      run = control$run,
      maxfitness = Inf,
      popSize = control$popSize
    )

    as.integer(result@solution[1,])
  }

  set_seriation_method("dist", "GA",
    GA_helper, "Genetic Algorithm")
}
