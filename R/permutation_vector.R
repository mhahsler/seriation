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



## S3 permutation and permutations classes
## permutations consists of instances of permutation

## permutation_vector

## constructor (NA is identity vector)
ser_permutation_vector <- function(x, method = NULL) {
  if(inherits(x, "ser_permutation_vector")) return(x)
  
  ## make sure it's an integer vector
  if(is.vector(x) && !is.integer(x)) x <- as.integer(x)

  if(.is_identity_permutation(x)) attr(x, "method") <- "identity permutation"
  if(!is.null(method)) attr(x, "method") <- method
    
  class(x) <- c("ser_permutation_vector", class(x))
  .valid_permutation_vector(x)
  x
}

## accessors
get_order <- function(x, ...) UseMethod("get_order")

get_order.ser_permutation_vector <- function(x, ...) NextMethod()

get_order.hclust <- function(x, ...) x$order

get_order.dendrogram <- function(x, ...) order.dendrogram(x)


get_order.integer <- function(x, ...) {
  o <- as.integer(x)
  if(.is_identity_permutation(x)) stop("Cannot get order from identity permutation.")
  #names(o) <- names(x)[o]
  o
}

## returns the order of objects (index of first, second, etc. object)
get_order.default <- function(x, ...) 
  stop(gettextf("No permutation accessor implemented for class '%s'. ",
    class(x)))

## returns for each object its rank (rank of first, second, etc. object)
get_rank <- function(x, ...) order(get_order(x, ...))

get_permutation_matrix <- function(x, ...) 
  permutation_vector2matrix(get_order(x, ...))

## c will create a ser_permutation!
c.ser_permutation_vector <- function(..., recursive = FALSE) 
  do.call("ser_permutation", list(...)) 

## convert to permuation matrix
permutation_vector2matrix <- function(x) {
  x <- get_order(x)
  .valid_permutation_vector(x)
  
  n <- length(x)
  pm <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n) pm[i, x[i]] <- 1
  pm
}

permutation_matrix2vector <- function(x) {
  .valid_permutation_matrix(x)
  o <- apply(x, MARGIN = 1, FUN = function(r) which(r==1))
  o
}

## reverse
rev.ser_permutation_vector <- function(x) {
  if(is(x, "hclust")) { x$order <- rev(x$order); x } 
  else ser_permutation_vector(rev(get_order(x)), method=get_method(x))
}


## currently method is an attribute of permutation
get_method <- function(x, printable = FALSE) {
  method <- attr(x, "method")
  
  if(printable && is.null(method)) method <- "unknown"
  method
}


## print et al
length.ser_permutation_vector <- function(x) 
  if(!.is_identity_permutation(x)) length(get_order(x)) else 0L 

print.ser_permutation_vector <-
  function(x, ...)
  {
    writeLines(
      c(gettextf("object of class %s",
        paste(sQuote(class(x)), collapse = ", ")),
        gettextf("contains a permutation vector of length %d", length(x)),
        gettextf("used seriation method: '%s'",
          get_method(x, printable = TRUE))))
    invisible(x)
  }

## fake summary (we dont really provide a summary, 
## but summary produces now a reasonable result --- same as print)
summary.ser_permutation_vector <- function(object, ...) {
  object
}


## helpers
.is_identity_permutation <- function (x) 
  if(is(x, "integer") && length(as.vector(x))==1 && is.na(x[1])) TRUE else FALSE

.valid_permutation_vector <- function(x) {
  ## identity vector
  if(.is_identity_permutation(x)) return()
  
  perm <- get_order(x)
  valid <- TRUE
  
  tab <- table(perm)
  if(any(tab != 1)) valid <- FALSE
  if(length(tab) != length(perm) 
    || any(names(tab) != sequence(length(perm)))) valid <- FALSE
  
  if(!valid) stop("Invalid permutation vector!\nVector: ", 
    paste(perm, collapse=", "))
}

.valid_permutation_matrix <- function(x) {
  if(any(rowSums(x)!=1) || any(colSums(x)!=1) || any(x!=1 & x!=0)) 
    stop("Not a valid permutation matrix")
}
