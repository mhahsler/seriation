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
function(x, order = NULL, method = NULL, force_loss = FALSE, ...)
  .criterion_array_helper(x, order, method, "matrix", force_loss)

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

### A MEASURE OF EFFECTIVENESS FOR THE MOMENT ORDERING ALGORITHM
### by Deutsch & Martin (1971)
### Correlation coefficient R for matrices.
criterion_R_matrix <- function(x, order, ...) {
  if(!is.null(order)) x <- permute(x, order)

  M <- nrow(x)
  N <- ncol(x)

  ## total sum
  T <- sum(x)

  ## X_i = i/M; Y_j = j/N
  X_i <- (1:M)/M
  Y_j <- (1:N)/N

  ## X_bar = 1/T sum_i,j a_ij X_i
  X_bar <- 1/T * sum(crossprod(x, X_i))
  ## Y_bar = 1/T sum_i,j a_ij Y_j
  Y_bar <- 1/T * sum(crossprod(t(x), Y_j))

  ## S_X2 = 1/(T-1) sum_i,j a_ij (X_i - X_bar)^2
  S_X2 <- 1/(T-1) * sum(crossprod(x, (X_i-X_bar)^2))
  ## S_Y2 = 1/(T-1) sum_i,j a_ij (Y_j - Y_bar)^2
  S_Y2 <- 1/(T-1) * sum(crossprod(t(x), (Y_j-Y_bar)^2))

  ## S_XY = 1/(T-1) sum_i,j a_ij  (X_i - X_bar) (Y_j - Y_bar)
  S_XY <- 1/(T-1) * sum(x * outer(X_i-X_bar, Y_j-Y_bar))

  ## R = S_XY/(S_X S_Y)
  S_XY/(sqrt(S_X2)*sqrt(S_Y2))
}

## register built-ins
set_criterion_method("matrix", "ME", criterion_ME,
    "Measure of effectiveness", TRUE)
set_criterion_method("matrix", "Cor_R", criterion_R_matrix,
    "Correlation Coefficient R", TRUE)

set_criterion_method("matrix", "Moore_stress", criterion_stress_moore,
    "Stress (Moore neighborhood)", FALSE)
set_criterion_method("matrix", "Neumann_stress", criterion_stress_neumann,
    "Stress (Neumann neighborhood)", FALSE)
