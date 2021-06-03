#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
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


## Seriation generic and default method.
seriate <- function(x, ...)
  UseMethod("seriate")
seriate.default <- function(x, ...)
  stop(gettextf("seriate not implemented for class '%s'.",
    class(x)))

## Seriation methods db.
get_seriation_method <- function(kind, name) {
  if (missing(kind))
    method <- registry_seriate$get_entry(name = name)
  else
    method <- registry_seriate$get_entry(kind = kind, name = name)

  if (is.null(method))
    stop(
      "Unknown seriation method ", 
      name,
      " for data type ",
      kind,
      ". Check list_seriation_methods(\"",
      kind,
      "\")"
    )

  method
}

set_seriation_method <- function(kind,
  name,
  definition,
  description = NULL,
  control = list(),
  ...) {
  ## check formals
  if (!identical(names(formals(definition)),
    c("x", "control")))
    stop("Seriation methods must have formals 'x' and 'control'.")

  ## check if entry already exists
  r <- registry_seriate$get_entry(kind = kind, name = name)
  if (!is.null(r) && r$name == name) {
    warning(
      "Entry with name \"",
      name,
      "\" for kind \"",
      kind,
      "\" already exists! Modifying entry."
    )
    registry_seriate$modify_entry(
      kind = kind,
      name = name,
      fun = definition,
      description = description,
      control = control
    )
  } else {
    registry_seriate$set_entry(
      name = name,
      kind = kind,
      fun = definition,
      description = description,
      control = control
    )
  }
}

list_seriation_methods <- function(kind) {
  if (missing(kind)) {
    kinds <- unique(sort(as.vector(
      sapply(registry_seriate$get_entries(), "[[", "kind")
    )))

    sapply(
      kinds,
      FUN = function(k)
        list_seriation_methods(k)
    )

  } else{
    sort(as.vector(sapply(
      registry_seriate$get_entries(kind = kind), "[[", "name"
    )))
  }
}

### deprecated
show_seriation_methods <- function(kind) {
  warning("Function is deprecated use: get_seriation_method() instead!")
  if (missing(kind))
    m <- registry_seriate$get_entries()
  else
    m <- registry_seriate$get_entries(kind = kind)
  m[sort(names(m))]
}
