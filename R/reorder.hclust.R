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

reorder.hclust <- function(x, dist, method = "OLO", ...) {
  method <- match.arg(tolower(method), choices = c("olo", "gw"))
  
  ## no reordering for less than 3 objects!
  if(length(x$order)<3) return(x)
  
  switch(method, 
      olo = .seriate_optimal(x, dist),
      gw  = .seriate_gruvaeus(x, dist)
    )
}

## wrapper for reorder.hclust in gclus
.seriate_gruvaeus <- function(hclust, dist)
    gclus::reorder.hclust(hclust, dist)

## wrapper to the optimal leaf ordering algorithm
##
## ceeboo 2005
.seriate_optimal <- function(hclust, dist) {
    ## check hclust
    merge <- hclust$merge
    if (!is.matrix(merge))
        stop("Component 'merge' of argument 'hclust' must be a matrix.")
    if (length(dim(merge)) != 2)
        stop("Component 'merge' of argument 'hclust' is invalid.")
    if (dim(merge)[1] != attr(dist,"Size")-1)
        stop("Argument 'dist' and component 'merge' of argument 'hclust' do not conform.")
    mode(merge) <- "integer"
    
    obj <- .Call("order_optimal", dist, merge)
    
    names(obj) <- c("merge","order","length")
    ##names(obj$order) <- attr(dist,"Labels")
    hclust$merge <- obj$merge
    hclust$order <- obj$order

    hclust
}
