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


#' @include seriate_TSP.R
.bea_tsp_contr <- .tsp_control

seriate_matrix_bea_tsp <-
  function(x, control, margin = seq_along(dim(x))) {
    if (any(x < 0))
      stop("Requires a nonnegative matrix.")

    # single objects do not work so we skip them
    if (1L %in% margin) {
      if (nrow(x) > 1L) {
        criterion <- as.dist(tcrossprod(x))
        row <- seriate(max(criterion) - criterion,
                       method = "TSP",
                       control = control)[[1]]
      } else {
        row <- 1L
      }
      attr(row, "method") <- "BEA_TSP"
    } else
      row <- NA

    if (2L %in% margin) {
      if (ncol(x) > 1L) {
        criterion <- as.dist(crossprod(x))
        col <- seriate(max(criterion) - criterion,
                       method = "TSP",
                       control = control)[[1]]
      } else {
        col <- 1L
      }
      attr(col, "method") <- "BEA_TSP"
    } else
      col <- NA

    list(row = row, col = col)
  }


seriate_matrix_bea <- function(x, control = NULL, margin = NULL) {
  control <- .get_parameters(control, list())

  ### BEA is just cheapest insertion
  control <- list(method = "cheapest_insertion",
                  two_opt = FALSE, rep = 1,
                  verbose = control$verbose)

  seriate_matrix_bea_tsp(x, control = control, margin = margin)
}

## register methods
set_seriation_method(
  "matrix",
  "BEA",
  seriate_matrix_bea,
  "Bond Energy Algorithm (BEA; McCormick 1972) to maximize the Measure of Effectiveness of a non-negative matrix.",
  list(),
  optimizes = .opt("ME", "Measure of effectiveness"),
  randomized = TRUE
)

set_seriation_method(
  "matrix",
  "BEA_TSP",
  seriate_matrix_bea_tsp,
  "Use a TSP to optimize the Measure of Effectiveness (Lenstra 1974).",
  .bea_tsp_contr,
  optimizes = .opt("ME", "Measure of effectiveness"),
  randomized = TRUE
)
