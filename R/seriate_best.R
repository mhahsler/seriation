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
#' Functions to automatically try different seriation methods or
#' rerun randomized methods several times to find the best order
#' given a criterion measure.
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
#' @param fast_only logical; exclude slow default seriation methods.
#' @param rep number of times to repeat the randomized seriation algorithm.
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
#' d_robinson <- as.dist(random.robinson(50, noise = .2, pre = TRUE))
#' m_iris <- as.matrix(iris[ sample(seq(nrow(iris))),-5])
#'
#' # find best seriation order
#' o <- seriate_best(d_supreme, verbose = TRUE)
#' pimage(d_supreme, o)
#'
#' o <- seriate_best(d_robinson, verbose = TRUE)
#' pimage(d_robinson, o)
#'
#' o <- seriate_best(m_iris, verbose = TRUE)
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
#' @export
seriate_best <- function(x,
                         methods = NULL,
                         criterion = NULL,
                         verbose = FALSE,
                         fast_only = TRUE,
                         ...) {
  ### data.frame/table?
  type <- class(x)[[1]]
  if (type %in% c("table", "data.frame"))
    type <- "matrix"

  # set some default methods
  if (is.null(methods)) {
    if (type == "dist") {
      methods <- c(
        #"ARSA", ## LS
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
      if (!fast_only)
        methods <- append("ARSA", methods)
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

  if (verbose)
    cat("Criterion:", criterion, "\n")
  cat("Performing: ")

  os <- sapply(
    methods,
    FUN = function(m) {
      if (verbose)
        cat(m, " ")
      tm <- system.time(o <- seriate(x, m, ...))
      attr(o, "time") <- tm[1] + tm[2]
      attr(o, "criterion") <- criterion(x, o, criterion,
                                        force_loss = TRUE)
      o
    },
    simplify = FALSE
  )

  if (verbose) {
    cat("\n\nResults:\n")
    df <- data.frame(
      method = names(os),
      criterion = sapply(os, attr, "criterion"),
      secs = sapply(os, attr, "time"),
      row.names = NULL
    )

    df <- df[order(df$criterion),]
    print(df)
    cat("\n")
  }

  os[[which.min(sapply(os, attr, "criterion"))]]
}

#' @rdname seriate_best
#' @export
seriate_rep <- function(x,
                        method,
                        criterion = NULL,
                        rep = 10L,
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
    warning("Specified seriation method is not randomized. Running it only once!")
    rep <- 1L
  }

  if (verbose)
    cat("Criterion:", criterion, "\nPerforming", rep, "tries\n")

  control <- list(...)
  r <- replicate(rep, seriate(x, method, control), simplify = FALSE)

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
      "Found orders with",
      criterion,
      "in the range:" ,
      min(cs),
      "to",
      max(cs),
      "\nReturning best."
    )


  o
}
