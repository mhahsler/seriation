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
#' It may not be know which seriation method produces the best result and
#' heuristics may produce unstable results.
#' `seriate_best()` and `seriate_rep()` automatically try different seriation methods or
#' rerun randomized methods several times to find the best and order
#' given a criterion measure.
#'
#' Rerun different seriation methods to find the best solution given a criterion
#' measure. Non-stochastic methods are automatically run only once.
#'
#' Support for parallel execution is provided using the [`foreach`] package. To
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
#' @param criterion a character string with the [criterion] to optimize.
#' @param verbose logical; show progress and results for different methods
#' @param rep number of times to repeat the randomized seriation algorithm.
#' @param parallel logical; perform replications in parallel.
#'      Uses `[foreach]` if a
#'      DoPar backend (e.g., `doParallel`) is rgistered.
#' @param ... further arguments are passed on (e.g., as `control`)
#'
#' @return Returns an object of class [ser_permutation].
#'
#' @author Michael Hahsler
#'
#' @keywords optimize cluster
#' @examples
#' # prepare some datasets (two distance matrices and a data matrix)
#' data(SupremeCourt)
#' d_supreme <- as.dist(SupremeCourt)
#' m_iris <- as.matrix(iris[ sample(seq(nrow(iris))),-5])
#'
#' # find best seriation order (tries by by default several fast methods)
#' o <- seriate_best(d_supreme, verbose = TRUE)
#' o
#' pimage(d_supreme, o)
#'
#' # specify the criterion
#' o <- seriate_best(d_supreme, criterion = "Path_length", verbose = TRUE)
#' o
#' pimage(d_supreme, o)
#'
#' # find best seriation for a matrix
#' o <- seriate_best(m_iris, verbose = TRUE)
#' o
#' pimage(m_iris, o, prop = FALSE)
#'
#' # run randomized algorithms several times. Repetition information
#' # is returned as attributes
#' o <- seriate_rep(d_supreme, "QAP_2SUM", rep = 10, verbose = TRUE)
#'
#' attr(o, "criterion")
#' hist(attr(o, "criterion_distribution"))
#' pimage(d_supreme, o)
#'
#' o <- seriate_rep(m_iris, "BEA_TSP", rep = 10, verbose = TRUE)
#'
#' attr(o, "criterion")
#' hist(attr(o, "criterion_distribution"))
#' pimage(m_iris, o, prop = FALSE)
#'
#' \dontrun{
#' # Using parallel execution
#'
#' library(doParallel)
#' registerDoParallel(cores = detectCores() - 1L)
#' o <- seriate_best(m_iris, verbose = TRUE, rep = 10)
#' o
#'
#'
#'
#' }
#' @export
seriate_best <- function(x,
                         methods = NULL,
                         criterion = NULL,
                         rep = 1L,
                         parallel = TRUE,
                         verbose = FALSE,
                         ...) {
  ### data.frame/table?
  type <- class(x)[[1]]
  if (type %in% c("table", "data.frame"))
    type <- "matrix"

  # set some default methods
  if (is.null(methods)) {
    if (type == "dist") {
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
    else if (type == "matrix")
      methods <- c("BEA_TSP", "PCA", "Heatmap", "PCA_angle")
    else
      stop("Currently only seriation for dist and matrix are supported.")
  }
  if (is.null(criterion)) {
    if (type == "dist")
      criterion <- "Gradient_weighted"
    else
      criterion <- "Moore_stress"
  }
  criterion <- get_criterion_method(type, criterion)$name

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
    df <- df[order(df$criterion), ]

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
                        method,
                        criterion = NULL,
                        rep = 10L,
                        parallel = TRUE,
                        verbose = TRUE,
                        ...) {
  type <- class(x)[[1]]
  if (type %in% c("table", "data.frame"))
    type <- "matrix"

  if (is.null(criterion)) {
    if (type == "dist")
      criterion <- "Gradient_weighted"
    else
      criterion <- "Moore_stress"
  }

  m <- get_seriation_method(type, method)
  method <- m$name
  if (!m$randomized && rep > 1L) {
    #message("Specified seriation method is not randomized. Running it only once!")
    rep <- 1L
  }

  if (verbose)
    cat("Criterion:", criterion, "\nTries", rep, " ")

  control <- list(...)
  #r <- replicate(rep, { if (verbose) cat("."); seriate(x, method, control) },
  #               simplify = FALSE)

  # r <- times(rep) %dopar% { list(seriate(x, method, control)) }

  dopar <-
    ifelse(foreach::getDoParRegistered() &&
             parallel && rep > 1L, `%dopar%`, `%do%`)

  r <-
    dopar(times(rep), {
      if (verbose)
        cat(".")
      list(seriate(x, method, control))
    })

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

  if (verbose)
    cat(
      "\nFound orders with",
      criterion,
      "in the range" ,
      min(cs),
      "to",
      max(cs),
      "- returning best\n"
    )

  o
}
