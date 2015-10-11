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



## Criterion for the quality of a permutation of a array

.criterion_array_helper <-
function(x, order = NULL, method = NULL, datatype = "array")
{
    ## check order
    if(!is.null(order)){
        if(!inherits(order, "ser_permutation")) 
            stop("Argument 'order' has to be of class 'ser_permutation'.")
        .check_matrix_perm(x, order)
    }

    ## get methods
    if(is.null(method)) method <- list_criterion_methods(datatype)
    method <- lapply(method, function(m) get_criterion_method(datatype, m))

    sapply(method,
        function(m) structure(m$fun(x, order), names=m$name))
}

criterion.array <-
function(x, order = NULL, method = NULL, ...)
    .criterion_array_helper(x, order, method, "array")

## methods

## register built-ins
