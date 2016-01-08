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



## seriate general arrays

.seriate_array_helper <- function(x, method = "PCA", control = NULL,
    margin = seq(length(dim(x))), datatype = "array", defmethod, ...){

    ## add ... to control
    control <- c(control, list(...))

    ## margin 1...rows, 2...cols, ...
    #if(is.null(method)) method <- "PCA"
    #else
    if(!is.character(method) || (length(method) != 1L))
      stop("Argument 'method' must be a character string.")

    method <- get_seriation_method(datatype, method)



    order <- method$fun(x, control)

    perm <- do.call("ser_permutation",
      unname(lapply(order, "ser_permutation_vector", method$name))
    )

    perm[margin]
  }

seriate.array <- function(x, method = "PCA", control = NULL,
    margin = seq(length(dim(x))), ...)
    .seriate_array_helper(x, method, control, margin,
      datatype = "array", defmethod = NA,...)
## we currently have no method and therefore also no default method!


## methods
## Identity is defined in seriate.matrix.R
## no other methods available right now

## register methods
## no methods available right now
