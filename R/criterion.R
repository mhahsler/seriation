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



## Criterion generic.
criterion <- function(x, order = NULL, method = NULL, ...)
  UseMethod("criterion")

## Criterion method registry.

## <NOTE>
## For criterion() methods, argument 'method' really allows selecting
## *several* methods ... should perhaps be called 'methods'?
## We thus have a getter which returns a named list of methods from the
## registry, and a setter for single methods.
## </NOTE>


set_criterion_method <- function(kind, name, fun, 
  description = NULL, merit = NA, ...) {
  ## check formals
  ##if(!identical(names(formals(definition)),
  ##              c("x", "order", "...")))
  ##    stop("Criterion methods must have formals 'x', 'order', and '...'.")
  
  registry_criterion$set_entry(
    kind = kind, name=name, fun = fun, 
    description = description, merit = merit)
}

get_criterion_method <- function(kind, name)
  registry_criterion$get_entry(kind=kind, name=name)

list_criterion_methods <- function(kind)
  as.vector(sapply(registry_criterion$get_entries(kind=kind), "[[", "name"))

show_criterion_methods <- function(kind)
  registry_criterion$get_entries(kind=kind)
