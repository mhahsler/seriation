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
criterion <- function(x, order = NULL, method = NULL, force_loss = FALSE, ...)
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

  ## check if criterion is already in registry
  r <- registry_criterion$get_entry(kind=kind, name=name)
  if(!is.null(r) && r$name==name) {
    warning("Entry with name ", name, " already exists! Modifying entry.")
    registry_criterion$modify_entry(kind=kind, name=name, fun=fun,
      description=description, merit = merit)
  } else {
    registry_criterion$set_entry(
      kind = kind, name=name, fun = fun,
      description = description, merit = merit)
  }

}

get_criterion_method <- function(kind, name) {
  method <- registry_criterion$get_entry(kind=kind, name=name)
  if(is.null(method))
    stop("Unknown criterion. Check list_criterion_methods(\"", kind, "\")")

  method
}

list_criterion_methods <- function(kind){
  if(missing(kind)) m <- registry_criterion$get_entries()
  else m <- registry_criterion$get_entries(kind=kind)

  sort(as.vector(sapply(m, "[[", "name")))
}

show_criterion_methods <- function(kind) {
  if(missing(kind)) m <- registry_criterion$get_entries()
  else m <- registry_criterion$get_entries(kind=kind)
  m[sort(names(m))]
}
