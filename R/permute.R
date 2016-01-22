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


## generic
permute <-
  function(x, order, ...)
    UseMethod("permute")

## methods
##permute.default <- function(x, order)
##stop(paste("\npermute not implemented for class: ", class(x)))
permute.default   <- function(x, order, ...) .permute_kd(x, order, ...)
permute.array     <- function(x, order, ...) .permute_kd(x, order, ...)
permute.matrix    <- function(x, order, ...) .permute_kd(x, order, ...)
permute.numeric   <- function(x, order, ...) .permute_1d(x, order, ...)
permute.character <- function(x, order, ...) .permute_1d(x, order, ...)
permute.list      <- function(x, order, ...) .permute_1d(x, order, ...)

## special cases
permute.dist <- function(x, order, ...){
  .nodots(...)

  if(!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  if(.is_identity_permutation(order[[1]])) return(x)

  .check_dist_perm(x, order)

  .rearrange_dist(x, get_order(order, 1))
}

permute.data.frame <- function(x, order, ...){
  .nodots(...)

  if(!inherits(order, "ser_permutation_vector"))
    order <- ser_permutation(order)

  if(length(order) != 1L)
    stop("dimensions do not match")

  perm <- get_order(order[[1L]])
  if(nrow(x) != length(perm))
    stop("some permutation vectors do not fit dimension of data")

  x[perm,]
}

permute.dendrogram <- function(x, order, ...) {
  .nodots(...)

  if(length(get_order(order)) != nobs(x))
    stop("Length of order and number of leaves in dendrogram do not agree!")

  x <- dendextend::rotate(x, order = match(get_order(order), order.dendrogram(x)))

  if(any(order.dendrogram(x) != get_order(order)))
    warning("Dendrogram cannot be perfectly reordered! Using best approximation.")

  x
}

permute.hclust <- function(x, order, ...) {
  nd <- as.hclust(permute(as.dendrogram(x), order, ...))
  x$merge <- nd$merge
  x$height <- nd$height
  x$order <- nd$order

  x
}

## helper
.check_dist_perm <- function(x, order){
  if(length(order) != 1L)
    stop("dimensions do not match")

  if(attr(x, "Size") != length(get_order(order, 1)))
    stop("some permutation vectors do not fit dimension of data")

  ## check dist
  if(attr(x, "Diag") || attr(x, "Upper"))
    stop("'dist' with diagonal or upper triangle matrix not implemented")
}

.check_matrix_perm <- function(x, order){
  if(length(dim(x)) != length(order))
    stop("dimensions do not match")
  if(any(dim(x) != sapply(order, length)))
    stop("some permutation vectors do not fit dimension of data")
}

.permute_kd <- function(x, order, ...){
  .nodots(...)

  if(!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  ## deal with identity permutations
  todo <- which(sapply(order, .is_identity_permutation))
  for(i in todo) order[[i]] <- ser_permutation_vector(seq(dim(x)[i]))

  .check_matrix_perm(x, order)

  perm <- lapply(order, get_order)
  do.call("[", c(list(x), perm, drop=FALSE))
}

.permute_1d <- function(x, order, ...) {
  .nodots(...)

  if(!inherits(order, "ser_permutation"))
    order <- ser_permutation(order)

  if(length(order) != 1)
    stop("dimensions do not match!")

  if(.is_identity_permutation(order[[1]])) return(x)

  perm <- get_order(order, 1)
  if(length(x) != length(perm))
    stop("some permutation vectors do not fit dimension of data!")

  x[perm]
}


## if we used proxy we would say:
#.rearrange_dist <- function (x, order) x[[order]]

.rearrange_dist <- function (x, order){
  ## make C call
  mode(x) <- "double"
  ## as.dist seems to make Size numeric and not integer!
  attr(x, "Size") <- as.integer(attr(x, "Size"))
  mode(order) <- "integer"

  d <- .Call("reorder_dist", x, order)

  labels <- if(is.null(labels(x)))
    NULL
  else
    labels(x)[order]

  structure(d,
    class   = "dist",
    Size    = length(order),
    Labels  = labels,
    Diag    = FALSE,
    Upper   = FALSE,
    method  = attr(x, "method")
  )
}
