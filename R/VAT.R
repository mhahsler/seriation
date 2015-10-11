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

## calculate path distance from iVAT using a modified version fo Floyd's alg.
## d_ij = smallest value of the largest values of all possible paths between i and j 
path_dist <- function(x) {
  
  #A <- as.matrix(x)
  #n <- nrow(A)
  #for(k in 1:n)
  # for(i in 1:n)
  #   for(j in 1:n) 
  #      if(max(A[i,k], A[k,j]) < A[i,j]) A[i,j] <- max(A[i,k], A[k,j])
  #d <- as.dist(A)
  
  ## make C call
  m <- as.matrix(x)
  
  if(any(is.na(m))) stop("NAs not allowed in x.")
  if(any(m<0)) stop("Negative values not allowed in x.")
  mode(m) <- "double"

  ## replace Inf with large number
  m[is.infinite(m)] <- .Machine$double.xmax
  
  if(any(m<0)) stop("Negative values not allowed in x.")
  
  m <- .Call("pathdist_floyd", m, PACKAGE = "seriation")
  as.dist(m)
}

VAT <- function(x, ...) {
  if(!is(x, "dist")) stop("x needs to be of class 'dist'!")
  pimage(x, seriate(x, "VAT"), ...)  
}

iVAT <- function(x, ...) {
  if(!is(x, "dist")) stop("x needs to be of class 'dist'!")
  x <- path_dist(x)
  pimage(x, seriate(x, "VAT"), ...)  
}