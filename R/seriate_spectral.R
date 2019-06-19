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


## Spectral Seriation
## Ding, C. and Xiaofeng He (2004): Linearized cluster assignment via
## spectral orderingProceedings of the Twenty-first.
## International Conference on Machine learning (ICML '04)

## Minimizes: sum_{i,j} (i-j)^2 * d_{pi_i,pi_j}

seriate_dist_spectral <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  ### calculate Laplacian
  W <- 1/(1+as.matrix(x))
  D <- diag(rowSums(W))
  L <- D - W

  ## The Fiedler vector is the eigenvector with the smallest eigenvalue
  ## eigen reports eigenvectors/values in decreasing order
  q <- eigen(L)
  fiedler <- q$vectors[,ncol(W)-1L]
  o <- order(fiedler)
  names(o) <- names(x)[o]
  o
}

seriate_dist_spectral_norm <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  ### calculate normalized Laplacian
  W <- 1/(1+as.matrix(x))
  D_sqrt<- diag(rowSums(1/W^.5))
  L <- D_sqrt %*% W %*% D_sqrt

  z <- eigen(L)$vectors
  q <- D_sqrt %*% z

  ## look for the vector with the largest eigenvalue
  largest_ev <- q[,2L]
  o <- order(largest_ev)
  names(o) <- names(x)[o]
  o
}


set_seriation_method("dist", "Spectral", seriate_dist_spectral,
  "Spectral seriation (Ding and He 2004)  uses a relaxation to minimize the 2-Sum Problem (Barnard, Pothen, and Simon 1993). It uses the order of the Fiedler vector of the similarity matrix's Laplacian.")
set_seriation_method("dist", "Spectral_norm", seriate_dist_spectral_norm,
  "Spectral seriation (Ding and He 2004)  uses a relaxation to minimize the 2-Sum Problem (Barnard, Pothen, and Simon 1993). It uses the order of the Fiedler vector of the similarity matrix's normalized Laplacian.")
