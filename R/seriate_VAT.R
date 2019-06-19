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



## VAT: a tool for visual assessment of (cluster) tendency
## Bezdek, J.C., Hathaway, R.J.
## Proceedings of the 2002 International Joint Conference on
## Neural Networks, 2002. IJCNN '02. (Volume:3)
seriate_dist_VAT <- function(x, control = NULL) {
  #param <- .get_parameters(control, NULL)
  .get_parameters(control, NULL)

  D <- as.matrix(x)
  N <- nrow(D)
  P <- rep(NA_integer_, N)
  I <- rep(FALSE, N)
  ### J is !I

  i <- which(D == max(D, na.rm = TRUE), arr.ind = TRUE)[1,1]
  P[1] <- i
  I[i] <- TRUE

  for(r in 2:N) {
    D2 <- D[I,!I, drop=FALSE]
    j <- which(D2 == min(D2, na.rm = TRUE), arr.ind = TRUE)[1,2]
    j <- which(!I)[j]
    P[r] <- j
    I[j] <- TRUE
  }

  names(P) <- labels(x)[P]
  P
}


set_seriation_method("dist", "VAT", seriate_dist_VAT,
  "Visual assesment of clustering tendency (Bezdek and Hathaway (2002). Creates an order based on Prim's algorithm for finding a minimum spanning tree (MST) in a weighted connected graph representing the distance matrix. The order is given by the order in which the nodes (objects) are added to the MST.")
