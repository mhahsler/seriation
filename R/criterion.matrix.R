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



## Criterion for the quality of a permutation of a matrix

criterion.matrix <-
function(x, order = NULL, method = NULL, ...)
  .criterion_array_helper(x, order, method, "matrix") 

## Bond energy (BEA)
criterion_ME <- function(x, order = NULL, ...) {
    ### ... unused

    if(any(x < 0))
        stop("Bond energy (ME) is only defined for nonnegative matrices.")
    
    n <- nrow(x)
    m <- ncol(x)
    
    if(!is.null(order)) x <- permute(x, order)
    
    mode(x) <- "single"

    energy <- .Fortran("energy",
                       n = n,
                       m = m,
                       b = x,
                       ener = as.single(0.0),
		       PACKAGE="seriation")
    
    0.5 * as.numeric(energy$ener)
}


## the interface to the stress functions allows for
## arbitrary subsetting (see the wrapper in C).
## (C) ceeboo 2005, 2006

.stress <- function(x, order, type="moore") {
    TYPE <- c(1,2)
    names(TYPE) <- c("moore", "neumann")
    if (inherits(x, "dist")) x <- as.matrix(x)
    if (!is.matrix(x))
        stop("Argument 'x' must be a matrix.")
    if (!is.double(x)) mode(x) <- "double"
    
    if(is.null(order)) {
        rows <- as.integer(1:dim(x)[1])
        cols <- as.integer(1:dim(x)[2])
    }else{
        rows <- get_order(order, 1)
        cols <- get_order(order, 2)
    }
    
    type <- as.integer(TYPE[type])
    
    x <- .Call("stress", x, rows, cols, type)
    
    ## does only half of the matrix!
    2 * x
}

criterion_stress_moore <- function(x, order, ...) .stress(x, order, "moore")
criterion_stress_neumann <- function(x, order, ...) .stress(x, order, "neumann")

## register built-ins
set_criterion_method("matrix", "ME", criterion_ME, 
    "Measure of effectiveness", TRUE)
set_criterion_method("matrix", "Moore_stress", criterion_stress_moore,
    "Stress (Moore neighborhood)", FALSE)
set_criterion_method("matrix", "Neumann_stress", criterion_stress_neumann,
    "Stress (Neumann neighborhood)", FALSE)
