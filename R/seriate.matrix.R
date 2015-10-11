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



## seriate matrices 

seriate.matrix <- function(x, method = "PCA", control = NULL, 
  margin = c(1,2), ...)
  .seriate_array_helper(x, method, control, margin, 
    datatype = "matrix", defmethod = "BEA_TSP", ...)

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
  
  if(any(x < 0)) stop("Requires a nonnegative matrix.")
  
  criterion <- as.dist(tcrossprod(x))
  row <- seriate(max(criterion)-criterion, 
    method = "TSP", control = control)[[1]]
  
  criterion <- as.dist(crossprod(x))
  col <- seriate(max(criterion)-criterion, 
    method = "TSP", control = control)[[1]]
  
  list(row = row, col = col)
}


## Bond Energy Algorithm (McCormick 1972)

seriate_matrix_bea <- function(x, control = NULL){
  control <- .get_parameters(control, list(
    istart = 0,
    jstart = 0,
    rep = 1
  ))
  
  if(any(x < 0)) stop("Requires a nonnegative matrix.")
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

## use the projection on the first pricipal component to determine the
## order
seriate_matrix_fpc <- function(x, control = NULL) {
  control <- .get_parameters(control, list(
    center = TRUE,
    scale. = FALSE,
    tol = NULL,
    verbose = FALSE
  ))
  
  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol
  verbose <- control$verbose
  
  pr <- prcomp(x, center = center, scale. = scale., tol = tol)
  scores <- pr$x[,1]
  row <- order(scores)
  if(verbose) cat("row: first principal component explains", 
    pr$sdev[1] / sum(pr$sdev)* 100,"%\n")
  
  pr <- prcomp(t(x), center = center, scale. = scale., tol = tol)
  scores <- pr$x[,1]
  col <- order(scores)
  if(verbose) cat("col: first principal component explains", 
    pr$sdev[1] / sum(pr$sdev)* 100,"%\n")
  
  names(row) <- rownames(x)[row]
  names(col) <- colnames(x)[col]
  
  list(row = row, col = col)
}

## Angle between the first 2 PCS. Fiendly (2002)
.order_angle <- function(x) {
  
  alpha <- atan2(x[,1], x[,2])
  o <- order(alpha)
  cut <- which.max(abs(diff(c(alpha[o], alpha[o[1]]+2*pi))))
  if(cut==length(o)) o
  else o[c((cut+1):length(o), 1:(cut))]
  
}


seriate_matrix_angle <- function(x, control = NULL) {
  control <- .get_parameters(control, list(
    center = TRUE,
    scale. = FALSE,
    tol = NULL
  ))
  
  center  <- control$center
  scale.  <- control$scale.
  tol     <- control$tol
  
  pr <- prcomp(x, center = center, scale. = scale., tol = tol)
  row <- .order_angle(pr$x[,1:2])
  
  pr <- prcomp(t(x), center = center, scale. = scale., tol = tol)
  col <- .order_angle(pr$x[,1:2])
  
  names(row) <- rownames(x)[row]
  names(col) <- colnames(x)[col]
  
  list(row = row, col = col)
}


seriate_matrix_identity <- function(x, control) {
  control <- .get_parameters(control, NULL)
  
  l <- lapply(dim(x), seq)
  for(i in 1:length(dim(x))) names(l[[i]]) <- labels(x)[[i]]
  l
}

seriate_matrix_random <- function(x, control) {
  control <- .get_parameters(control, NULL)
  
  l <- lapply(dim(x), FUN = function(l) sample(seq(l)))
  for(i in 1:length(dim(x))) names(l[[i]]) <- labels(x)[[i]][l[[i]]]
  l
}


## register methods
set_seriation_method("matrix", "BEA_TSP", seriate_matrix_bea_tsp,
  "TSP to maximize ME")
set_seriation_method("matrix", "BEA", seriate_matrix_bea,
  "Bond Energy Algorithm to maximize ME")
set_seriation_method("matrix", "PCA", seriate_matrix_fpc,
  "First principal component")
set_seriation_method("matrix", "PCA_angle", seriate_matrix_angle,
  "First two principal components (angle)")

set_seriation_method("matrix", "Identity", seriate_matrix_identity,
  "Identity permutation")
set_seriation_method("matrix", "Random", seriate_matrix_random,
  "Random permutation")

#set_seriation_method("matrix", "SVD", seriate_matrix_svd,
#    "Angles formed by first two eigenvectors")

## Array
set_seriation_method("array", "Identity", seriate_matrix_identity,
  "Identity permutation")
set_seriation_method("array", "Random", seriate_matrix_random,
  "Random permutation")

