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


## BEA FORTRAN code by Fionn Murtagh

bea <- function(a, istart=0, jstart=0)
{
    ## Permute rows and colums, using "bond energy algorithm".
    if(!is.matrix(a))
        stop("First input argument must be a matrix.\n")

    n <- nrow(a)
    m <- ncol(a)
    b    <- matrix(0.0, n, m)
    mode(a) <- "single"
    mode(b) <- "single"
    ib   <- integer(n)
    jb   <- integer(m)
    ifin <- integer(n)
    jfin <- integer(m)
    ener <- 0.0

    if (istart == 0) istart <- floor(runif(1,1,n))
    if (jstart == 0) jstart <- floor(runif(1,1,m))

    bea1  <- .Fortran("rbea",
                      n = as.integer(n),
                      m = as.integer(m),
                      a = as.matrix(a),             # input data
                      istart = as.integer(istart),  # 1st row placement
                      b = as.matrix(b),             # permuted array
                      ib = as.integer(ib),          # permuted order of rows
		      ifin = as.integer(ifin),      # for book-keeping
		      PACKAGE="seriation")
    a <- bea1$b
    
    bea2  <- .Fortran("cbea",
                      n = as.integer(n),
                      m = as.integer(m),
                      a = as.matrix(a),             # input data
                      jstart = as.integer(jstart),  # 1st col. placement
                      b = as.matrix(b),             # permuted array
                      jb = as.integer(jb),          # permuted order of cols.
		      jfin = as.integer(jfin),      # for book-keeping
		      PACKAGE="seriation")

    energ  <- .Fortran("energy",
                       n = as.integer(n),
                       m = as.integer(m),
                       b = as.matrix(bea2$b),
                       ener = as.single(ener),
		       PACKAGE="seriation")

    list(b = bea2$b, ib = bea1$ib, jb = bea2$jb, e = energ$ener)
}
