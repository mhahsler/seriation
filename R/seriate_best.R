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

#' Best Seriation
#'
#' Often the best seriation method for a particular dataset is not know and
#' heuristics may produce unstable results.
#' `seriate_best()` and `seriate_rep()` automatically try different seriation methods or
#' rerun randomized methods several times to find the best and order
#' given a criterion measure. `seriate_improve()` uses a local improvement strategy
#' to imporve an existing solution.
#'
#' `seriate_rep()` rerun a randomized seriation methods to find the best solution
#' given the criterion specified for the method in the registry.
#' A specific criterion can also be specified.
#' Non-stochastic methods are automatically only run once.
#'
#' `seriate_best()` runs a set of methods and returns the best result given a
#' criterion. Stochastic methods are automatically randomly restarted several times.
#'
#' `seriate_improve()` improves a seriation order using simulated annealing using
#' a specified criterion measure. It uses [seriate()] with method "`GSA`",
#' a reduced probability to accept bad moves, and a lower minimum temperature. Control
#' parameters for this method are accepted.
#'
#' **Criterion**
#'
#' If no criterion is specified, ten the criterion specified for the method in
#' the registry (see `[get_seriation_method()]`) is used. For methods with no
#' criterion in the registry (marked as "other"), a default method is used.
#' The defaults are:
#'
#' * `dist`: `"AR_deviations"` - the study in Hahsler (2007) has shown that this
#'  criterion has high similarity with most other criteria.
#'  * `matrix`: "Moore_stress"
#'
#' **Parallel Execution**
#'
#' Some methods support for parallel execution is provided using the [`foreach`] package. To
#' use parallel execution, a suitable backend needs to be registered (eee
#' the Examples section for using the `doParallel` package).
#'
#' @family seriation
#'
#' @param x the data.
#' @param method a character string with the name of the seriation method
#' (default: varies by data type).
#' @param methods a vector of character string with the name of the seriation
#' methods to try.
#' @param control a list of control options passed on to [seriate()].
#'      For `seriate_best()` control needs to be a named list of control lists
#'      with the names matching the seriation methods.
#' @param criterion `seriate_rep()` chooses the criterion specified for the
#'    method in the registry. A character string with the [criterion] to optimize
#'    can be specified.
#' @param verbose logical; show progress and results for different methods
#' @param rep number of times to repeat the randomized seriation algorithm.
#' @param parallel logical; perform replications in parallel.
#'      Uses `[foreach]` if a
#'      DoPar backend (e.g., `doParallel`) is rgistered.
#' @param ... further arguments are passed on to the [seriate()].
#'
#' @return Returns an object of class [ser_permutation].
#'
#' @author Michael Hahsler
#'
#' @keywords optimize cluster
#' @references
#' Hahsler, M. (2017): An experimental comparison of seriation methods for
#' one-mode two-way data. \emph{European Journal of Operational Research,}
#' \bold{257}, 133--143.
#' \doi{10.1016/j.ejor.2016.08.066}
#'
#' @examples
#' data(SupremeCourt)
#' d_supreme <- as.dist(SupremeCourt)
#'
#' # find best seriation order (tries by by default several fast methods)
#' o <- seriate_best(d_supreme, criterion = "AR_events")
#' o
#' pimage(d_supreme, o)
#'
#' # run a randomized algorithms several times. It automatically chooses the
#' # LS criterion. Repetition information is returned as attributes
#' o <- seriate_rep(d_supreme, "QAP_LS", rep = 5)
#'
#' attr(o, "criterion")
#' hist(attr(o, "criterion_distribution"))
#' pimage(d_supreme, o)
#'
#' \dontrun{
#' # Using parallel execution on a larger dataset
#' data(iris)
#' m_iris <- as.matrix(iris[sample(seq(nrow(iris))),-5])
#' d_iris <- dist(m_iris)
#'
#' library(doParallel)
#' registerDoParallel(cores = detectCores() - 1L)
#'
#' # seriate rows of the iris data set
#' o <- seriate_best(d_iris, criterion = "LS")
#' o
#'
#' pimage(d_iris, o)
#'
#' # improve the order to minimize RGAR instead of LS
#' o_improved <- seriate_improve(d_iris, o, criterion = "RGAR")
#' pimage(d_iris, o_improved)
#'
#' # available control parameters for seriate_improve()
#' get_seriation_method(name = "GSA")
#' }
#' @export
seriate_best <- function(x,
                         methods = NULL,
                         control = NULL,
                         criterion = NULL,
                         rep = 10L,
                         parallel = TRUE,
                         verbose = TRUE,
                         ...) {
  ### data.frame/table?
  kind <- get_seriation_kind(x)

  # set some default methods
  if (is.null(methods)) {
    if (kind == "dist") {
      methods <- c(
        "spectral",
        ## 2-Sum
        "MDS",
        ## Moore stress
        "QAP_2SUM",
        "QAP_BAR",
        "QAP_LS",
        "QAP_Inertia",
        "TSP",
        ## path length
        "OLO_average" ## restricted path length
      )
    }
    else if (kind == "matrix")
      methods <- c("BEA_TSP", "PCA", "Heatmap", "PCA_angle")
    else
      stop("Currently only seriation for dist and matrix are supported.")
  }

  if (is.null(criterion))
      criterion <- get_default_criterion(x)
  criterion <- get_criterion_method(kind, criterion)$name

  if (verbose) {
    cat("Criterion:", criterion, "\n")
    cat("Performing: ")
  }

  os <- sapply(
    methods,
    FUN = function(m) {
      if (verbose) {
        cat("\n")
        cat(m, " - ")
      }
      #tm <- system.time(o <- seriate(x, m, ...))
      tm <-
        system.time(
          o <-
            seriate_rep(
              x,
              m,
              control = control[[m]],
              verbose = verbose,
              criterion = criterion,
              rep = rep,
              parallel = parallel,
              ...
            )
        )
      attr(o, "time") <- tm[1] + tm[2]
      attr(o, "criterion") <- criterion(x, o, criterion,
                                        force_loss = TRUE)
      o
    },
    simplify = FALSE
  )

  if (verbose) {
    df <- data.frame(
      method = names(os),
      criterion = sapply(os, attr, "criterion"),
      secs = sapply(os, attr, "time"),
      row.names = NULL
    )
    df <- df[order(df$criterion),]

    cat("\nResults (first was chosen):\n")
    print(df)
    cat("\n")
  }

  os[[which.min(sapply(os, attr, "criterion"))]]
}

#' @rdname seriate_best
#' @importFrom foreach times `%dopar%` `%do%`
#' @export
seriate_rep <- function(x,
                        method = NULL,
                        control = NULL,
                        criterion = NULL,
                        rep = 10L,
                        parallel = TRUE,
                        verbose = TRUE,
                        ...) {
  if (is.null(method))
    method <- get_default_method(x)

  m <- get_seriation_method(get_seriation_kind(x), method)
  method <- m$name

  if (is.null(criterion))
    criterion <- m$optimizes

  if (is.na(criterion))
    criterion <- get_default_criterion(x)


  if (!m$randomized && rep > 1L) {
    rep <- 1L
    if (verbose)
      cat("Method not randomized. Running once")
  }

  if (verbose && rep > 1L) {
    cat("Tries", rep, " ")
  }

  #r <- replicate(rep, { if (verbose) cat("."); seriate(x, method, control) },
  #               simplify = FALSE)

  # r <- times(rep) %dopar% { list(seriate(x, method, control)) }

  dopar <-
    ifelse(foreach::getDoParRegistered() &&
             parallel && rep > 1L,
           `%dopar%`,
           `%do%`)

  r <-
    dopar(times(rep), {
      if (verbose)
        cat(".")
      list(seriate(x, method, control, ...))
    })

  if (verbose)
    cat("\n")

  cs <- sapply(
    r,
    FUN = function(o)
      criterion(x, o, criterion,
                force_loss = TRUE)
  )

  o <- r[[which.min(cs)]]
  attr(o, "criterion") <- min(cs)
  attr(o, "criterion_method") <- criterion
  attr(o, "criterion_distribution") <- as.vector(cs)


  if (verbose && rep > 1L)
    cat(
      "Found orders with",
      sQuote(criterion),
      "in the range" ,
      min(cs),
      "to",
      max(cs),
      "- returning best\n"
    )

  o
}

#' @rdname seriate_best
#' @param order a `ser_permutation` object for `x` or the name of a seriation method to start with.
#' @export
seriate_improve <- function(x,
                            order,
                            criterion = NULL,
                            control = NULL,
                            verbose = TRUE,
                            ...) {
  if (is.null(criterion))
    criterion <- get_default_criterion(x)

  criterion <-
    get_criterion_method(get_seriation_kind(x), criterion)$name

  if (is.null(control))
    control <- list()

  if (is.null(control$p_initial))
    control$p_initial <- 0.01 * 1e-6
  if (is.null(control$t_min))
    control$t_min <- 1e-12
  control$warmstart <- order
  control$criterion <- criterion
  control$verbose <- verbose

  seriate(x, "GSA", control = control, ...)
}
