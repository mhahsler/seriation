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


#' Registry for Seriation Methods
#'
#' A registry to manage methods used by [seriate()].
#'
#' The functions below are convenience function for the registry
#' \code{registry_seriate}.
#'
#' \code{list_seriation_method()} lists all available methods for a given data
#' type (\code{kind}). The result is a vector of character strings with the
#' short names of the methods. If \code{kind} is missing, then a list of
#' methods is returned.
#'
#' \code{get_seriation_method()} returns information (including the
#' implementing function) about a given method in form of an object of class
#' \code{"seriation_method"}.
#'
#' With \code{set_seriation_method()} new seriation methods can be added by the
#' user. The implementing function (\code{definition}) needs to have the formal
#' arguments \code{x, control}, where \code{x} is the data object and
#' \code{control} contains a list with additional information for the method
#' passed on from \code{seriate()}.  The implementation has to return a list of
#' objects which can be coerced into \code{ser_permutation_vector} objects
#' (e.g., integer vectors). The elements in the list have to be in
#' corresponding order to the dimensions of \code{x}.
#'
#' @import registry
#' @param kind the data type the method works on. For example, \code{"dist"},
#' \code{"matrix"} or \code{"array"}. If missing, then methods for any type are
#' shown.
#' @param name a short name for the method used to refer to the method in
#' [seriate()].
#' @param definition a function containing the method's code.
#' @param description a description of the method. For example, a long name.
#' @param control a list with control arguments and default values.
#' @param ... further information that is stored for the method in the
#' registry.
#' @returns
#' \code{list_seriation_method()} result is a vector of character strings with the
#' short names of the methods.
#'
#' \code{get_seriation_method()} returns a given method in form of an object of class
#' \code{"seriation_method"}.
#'
#' @author Michael Hahsler
#' @keywords misc
#' @examples
#'
#' # Registry
#' registry_seriate
#'
#' # List all seriation methods by type
#' list_seriation_methods()
#'
#' # List methods for matrix seriation
#' list_seriation_methods("matrix")
#'
#' get_seriation_method(name = "BEA")
#'
#' # Example for defining a new seriation method (reverse identity function for matrix)
#'
#' # 1. Create the seriation method
#' seriation_method_reverse <- function(x, control) {
#'    # return a list of order vectors, one for each dimension
#'    list(seq(nrow(x), 1), seq(ncol(x), 1))
#' }
#'
#' # 2. Register new method
#' set_seriation_method("matrix", "Reverse", seriation_method_reverse,
#'     description = "Reverse identity order", control = list())
#'
#' list_seriation_methods("matrix")
#' get_seriation_method("matrix", "reverse")
#'
#' # 3. Use the new seriation methods
#' seriate(matrix(1:12, ncol=3), "reverse")
#' @export
registry_seriate <- registry(registry_class = "seriation_registry",
  entry_class = "seriation_method")

registry_seriate$set_field("kind",
  type = "character",
  is_key = TRUE,
  index_FUN = match_partial_ignorecase)
registry_seriate$set_field("name",
  type = "character",
  is_key = TRUE,
  index_FUN = match_partial_ignorecase)
registry_seriate$set_field("fun", type = "function",
  is_key = FALSE)
registry_seriate$set_field("description", type = "character",
  is_key = FALSE)
registry_seriate$set_field("control", type = "list",
  is_key = FALSE)


print.seriation_method <- function(x, ...) {
  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    gettextf("description: %s", x$description)
  ))

  if (length(x$control) > 0) {
    writeLines("control (default values):")

    contr <- lapply(
      x$control,
      FUN =
        function(p)
          utils::capture.output(dput(p, control = list()))[1]
    )

    print(as.data.frame(contr))
  } else
    writeLines("control: no parameters registered.")

  invisible(x)
}

#' @rdname registry_seriate
#' @export
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

#' @rdname registry_seriate
#' @export
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
    # warning(
    #   "Entry with name \"",
    #   name,
    #   "\" for kind \"",
    #   kind,
    #   "\" already exists! Modifying entry."
    # )
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

#' @rdname registry_seriate
#' @export
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
#' @rdname registry_seriate
#' @export
show_seriation_methods <- function(kind) {
  warning("Function is deprecated use: get_seriation_method() instead!")
  if (missing(kind))
    m <- registry_seriate$get_entries()
  else
    m <- registry_seriate$get_entries(kind = kind)
  m[sort(names(m))]
}


#' Registry for Criterion Methods
#'
#' A registry to manage methods used by [criterion()] to calculate a criterion value given data and a
#' permutation.
#'
#' All methods below are convenience methods for the registry named
#' `registry_criterion`.
#'
#' `list_criterion_method()` lists all available methods for a given data
#' type (`kind`). The result is a vector of character strings with the
#' short names of the methods. If `kind` is missing, then a list of
#' methods is returned.
#'
#' `get_criterion_method()` returns information (including the
#' implementing function) about a given method in form of an object of class
#' `"criterion_method"`.
#'
#' With `set_criterion_method()` new criterion methods can be added by the
#' user. The implementing function (`fun`) needs to have the formal
#' arguments `x, order, ...`, where `x` is the data object, order is
#' an object of class [ser_permutation_vector] and `...` can contain
#' additional information for the method passed on from [criterion()]. The
#' implementation has to return the criterion value as a scalar.
#'
#' @name registry_criterion
#' @aliases registry_criterion registry
#' @param kind the data type the method works on. For example, `"dist"`,
#' `"matrix"` or `"array"`.
#' @param name a short name for the method used to refer to the method in the
#' function [criterion()].
#' @param fun a function containing the method's code.
#' @param description a description of the method. For example, a long name.
#' @param merit a boolean indicating if the criterion measure is a merit
#' (`TRUE`) or a loss (`FALSE`) measure.
#' @param ... further information that is stored for the method in the
#' registry.
#' @returns
#' `list_criterion_method()` results is a vector of character strings with the
#' short names of the methods.
#'
#' `get_criterion_method()` returns a given method in form of an object of class
#' `"criterion_method"`.
#' @author Michael Hahsler
#' @keywords misc
#' @examples
#' ## the registry
#' registry_criterion
#'
#' # List all criterion calculation methods by type
#' list_criterion_methods()
#'
#' # List methods for matrix
#' list_criterion_methods("matrix")
#'
#' get_criterion_method("dist", "AR_d")
#'
#' # Define a new method (sum of the diagonal elements)
#'
#' ## 1. implement a function to calculate the measure
#' criterion_method_matrix_foo <- function(x, order, ...) {
#' if(!is.null(order)) x <- permute(x,order)
#'     sum(diag(x))
#' }
#'
#' ## 2. Register new method
#' set_criterion_method("matrix", "DiagSum", criterion_method_matrix_foo,
#'     description = "Calculated the sum of all diagonal entries", merit = FALSE)
#'
#' list_criterion_methods("matrix")
#' get_criterion_method("matrix", "DiagSum")
#'
#' ## 3. use all criterion methods (including the new one)
#' criterion(matrix(1:9, ncol = 3))
#' @export
registry_criterion <-
  registry(registry_class = "criterion_registry",
    entry_class = "criterion_method")

registry_criterion$set_field("kind",
  type = "character",
  is_key = TRUE,
  index_FUN = match_partial_ignorecase)
registry_criterion$set_field("name",
  type = "character",
  is_key = TRUE,
  index_FUN = match_partial_ignorecase)
registry_criterion$set_field("fun", type = "function",
  is_key = FALSE)
registry_criterion$set_field("description", type = "character",
  is_key = FALSE)
registry_criterion$set_field("merit", type = "logical",
  is_key = FALSE)


print.criterion_method <- function(x, ...) {
  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    gettextf("description: %s", x$description),
    gettextf("merit:       %s", x$merit)
  ))
  invisible(x)
}

## <NOTE>
## For criterion() methods, argument 'method' really allows selecting
## *several* methods ... should perhaps be called 'methods'?
## We thus have a getter which returns a named list of methods from the
## registry, and a setter for single methods.
## </NOTE>

#' @rdname registry_criterion
#' @export
set_criterion_method <- function(kind,
  name,
  fun,
  description = NULL,
  merit = NA,
  ...) {
  ## check formals
  ##if(!identical(names(formals(definition)),
  ##              c("x", "order", "...")))
  ##    stop("Criterion methods must have formals 'x', 'order', and '...'.")

  ## check if criterion is already in registry
  r <- registry_criterion$get_entry(kind = kind, name = name)
  if (!is.null(r) && r$name == name) {
    warning("Entry with name ", name, " already exists! Modifying entry.")
    registry_criterion$modify_entry(
      kind = kind,
      name = name,
      fun = fun,
      description = description,
      merit = merit
    )
  } else {
    registry_criterion$set_entry(
      kind = kind,
      name = name,
      fun = fun,
      description = description,
      merit = merit
    )
  }

}

#' @rdname registry_criterion
#' @export
get_criterion_method <- function(kind, name) {
  method <- registry_criterion$get_entry(kind = kind, name = name)
  if (is.null(method))
    stop("Unknown criterion. Check list_criterion_methods(\"",
      kind,
      "\")")

  method
}

#' @rdname registry_criterion
#' @export
list_criterion_methods <- function(kind) {
  if (missing(kind)) {
    kinds <- unique(sort(as.vector(
      sapply(registry_criterion$get_entries(), "[[", "kind")
    )))

    sapply(
      kinds,
      FUN = function(k)
        list_criterion_methods(k)
    )

  } else{
    sort(as.vector(sapply(
      registry_criterion$get_entries(kind = kind), "[[", "name"
    )))
  }
}

### deprecated
#' @rdname registry_criterion
#' @export
show_criterion_methods <- function(kind) {
  warning("Function is deprecated use: get_criterion_method() instead!")
  if (missing(kind))
    m <- registry_criterion$get_entries()
  else
    m <- registry_criterion$get_entries(kind = kind)
  m[sort(names(m))]
}
