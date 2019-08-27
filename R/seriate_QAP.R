#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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

## QAP 2SUM seriation
seriate_dist_2SUM <- function(x, control = NULL) {
  ## param are passed on to QAP

  do.call(qap::qap, c(list(A = .A_2SUM(attr(x, "Size")),
    B = 1/(1+as.matrix(x))), control))
}

## QAP Linear seriation
seriate_dist_LS <- function(x, control = NULL) {
  ## param are passed on to QAP

  do.call(qap::qap, c(list(A = .A_LS(attr(x, "Size")),
    B = as.matrix(x)), control))
}

## QAP Inertia
seriate_dist_Inertia <- function(x, control = NULL) {
  ## param are passed on to QAP
  n <- attr(x, "Size")

  ## inertia uses the same A matrix as 2SUM
  ## we use n^2 since A needs to be positive
  do.call(qap::qap, c(list(A = n^2-.A_2SUM(n),
    B = as.matrix(x)), control))
}

## QAP BAR
.qap_bar_contr <- list(b = function(n) max(1, floor(n/5)))
seriate_dist_BAR <- function(x, control = NULL) {
  ## param are passed on to QAP

  if(is.null(control)) control <- .qap_bar_contr
  if(is.null(control$b)) control$b <- .qap_bar_contr$b

  .A_BAR <- function(n, b) {
    b <- floor(b)
    if(b<1 || b>=n) stop("b: needs to be 1<=b<n!")
    A <- b+1- outer(1:n,1:n, FUN = function(i,j) abs(i-j))
    A[A<0] <- 0
    A
  }

  n <- attr(x, "Size")
  if(is.function(control$b)) b <- control$b(n)
  else b <- floor(control$b)

  if(b<1 || b>n) stop("BAR bandwidth is not between 1 and n!")

  control$b <- NULL

  ## inertia uses the same A matrix as 2SUM
  do.call(qap::qap, c(list(A = .A_BAR(n, b = b),
    B = as.matrix(x)), control))
}


set_seriation_method("dist", "QAP_2SUM", seriate_dist_2SUM,
  "Quadratic assignment problem formulation for seriation solved using a simulated annealing solver to minimize the 2-Sum Problem criterion (Barnard, Pothen, and Simon 1993). Control arguments are passed to qap in package qap.")
set_seriation_method("dist", "QAP_LS", seriate_dist_LS,
  "Quadratic assignment problem formulation for seriation solved using a simulated annealing solver to minimize the Linear Seriation Problem (LS) criterion (Hubert and Schultz 1976).Control arguments are passed to qap in package qap.")
set_seriation_method("dist", "QAP_BAR", seriate_dist_BAR,
  "Quadratic assignment problem formulation for seriation solved using a simulated annealing solver to minimize the banded anti-Robinson form (BAR). Control argument b is the BAR bandwidth (either a function of n or a number). The default for b is n/5. Additional control arguments are passed to qap in package qap.", .qap_bar_contr)
set_seriation_method("dist", "QAP_Inertia", seriate_dist_Inertia,
  "Quadratic assignment problem formulation for seriation solved using a simulated annealing solver to minimize the Inertia criterion. Control arguments are passed to qap in package qap.")
