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


## Algorithm B
##  F. Murtagh (1985). Multidimensional Cluster Algorithms. Lectures
##  in Computational Statistics, Physica Verlag, pp. 15.
#
# this is actually just the same as BEA
#
#.seriate_matrix_murtagh <- function(x, control) {
#
#    if(any(x < 0)) stop("Requires a nonnegative matrix.")
#
#    criterion <- as.dist(tcrossprod(x))
#    row <- hclust_greedy(-criterion)$order
#    criterion <- as.dist(crossprod(x))
#    col <- hclust_greedy(-criterion)$order
#
#    list(row = row, col = col)
#}

seriate_matrix_bea_tsp <- function(x, control) {
  if (any(x < 0))
    stop("Requires a nonnegative matrix.")

  criterion <- as.dist(tcrossprod(x))
  row <- seriate(max(criterion) - criterion,
    method = "TSP",
    control = control)[[1]]

  criterion <- as.dist(crossprod(x))
  col <- seriate(max(criterion) - criterion,
    method = "TSP",
    control = control)[[1]]

  attr(row, "method") <- "BEA_TSP"
  attr(col, "method") <- "BEA_TSP"

  list(row = row, col = col)
}


## Bond Energy Algorithm (McCormick 1972)
.bea_contr <- list(istart = 0,
  jstart = 0,
  rep = 1)

seriate_matrix_bea <- function(x, control = NULL) {
  control <- .get_parameters(control, .bea_contr)

  if (any(x < 0))
    stop("Requires a nonnegative matrix.")
  istart <- control$istart
  jstart <- control$jstart
  rep  <- control$rep

  res <- replicate(rep, bea(x, istart = istart, jstart = jstart),
    simplify = FALSE)

  best <- which.max(sapply(res, "[[", "e"))
  res <- res[[best]]

  row <- res$ib
  col <- res$jb

  names(row) <- rownames(x)[row]
  names(col) <- colnames(x)[col]

  list(row = row, col = col)
}

## register methods
set_seriation_method(
  "matrix",
  "BEA",
  seriate_matrix_bea,
  "Bond Energy Algorithm (BEA; McCormick 1972) to maximize the Measure of Effectiveness of a non-negative matrix.",
  .bea_contr
)
set_seriation_method(
  "matrix",
  "BEA_TSP",
  seriate_matrix_bea_tsp,
  "Use a TSP to optimize the Measure of Effectiveness (Lenstra 1974). Control is passed on to the seriation method TSP."
)
