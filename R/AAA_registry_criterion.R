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
#' @name registry_for_criterion_methods
#' @family criterion
#'
#' @param kind the data type the method works on. For example, `"dist"`,
#' `"matrix"` or `"array"`.
#' @param name the name for the method used to refer to the method in the
#' function [criterion()].
#' @param names_only logical; return only the method name. `FALSE` returns
#'    also the method descriptions.
#' @param fun a function containing the method's code.
#' @param description a description of the method. For example, a long name.
#' @param merit a boolean indicating if the criterion measure is a merit
#' (`TRUE`) or a loss (`FALSE`) measure.
#' @param x an object of class "criterion_method" to be printed.
#' @param verbose logical; print a message when a new method is registered.
#' @param ... further information that is stored for the method in the
#' registry.
#' @returns
#' - `list_criterion_method()` results is a vector of character strings with the
#'   names of the methods used for `criterion()`.
#' - `get_criterion_method()` returns a given method in form of an object of class
#'   `"criterion_method"`.
#' @author Michael Hahsler
#' @seealso This registry uses [registry()] in package \pkg{registry}.
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
#' # get more description
#' list_criterion_methods("matrix", names_only = FALSE)
#'
#' # get a specific method
#' get_criterion_method(kind = "dist", name = "AR_d")
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

#' @rdname registry_for_criterion_methods
#' @export
list_criterion_methods <- function(kind, names_only = TRUE) {
  if (missing(kind)) {
    kinds <- unique(sort(as.vector(
      sapply(registry_criterion$get_entries(), "[[", "kind")
    )))

    sapply(
      kinds,
      FUN = function(k)
        list_criterion_methods(k, names_only = names_only)
    )

  } else{
    if (names_only)

      sort(as.vector(sapply(
        registry_criterion$get_entries(kind = kind), "[[", "name"
      )))
    else {
      l <- registry_criterion$get_entries(kind = kind)
      l[order(names(l))]
    }
  }
}


#' @rdname registry_for_criterion_methods
#' @export
get_criterion_method <- function(kind, name) {
  if (missing(kind))
    method <- registry_criterion$get_entry(name = name)
  else
    method <- registry_criterion$get_entry(kind = kind, name = name)

  if (is.null(method))
    stop("Unknown criterion. Check list_criterion_methods()")

  method
}

## <NOTE>
## For criterion() methods, argument 'method' really allows selecting
## *several* methods ... should perhaps be called 'methods'?
## We thus have a getter which returns a named list of methods from the
## registry, and a setter for single methods.
## </NOTE>

#' @rdname registry_for_criterion_methods
#' @export
set_criterion_method <- function(kind,
                                 name,
                                 fun,
                                 description = NULL,
                                 merit = NA,
                                 verbose = FALSE,
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

  if (verbose)
    message("Registering new seriation criteron ",
        sQuote(name),
        " for ",
        sQuote(kind))

}

#' @rdname registry_for_criterion_methods
#' @export
print.criterion_method <- function(x, ...) {
  extra_param <- setdiff(names(as.list(args(x$fun))), c("x", "order", "...", ""))

  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    strwrap(
      gettextf("description: %s", x$description),
      prefix = "             ",
      initial = ""
    ),
    gettextf("merit:       %s", x$merit)
  ))

  if (length(extra_param) > 0L)
    cat("parameters: ", paste(extra_param, collapse = ", "), "\n")

  invisible(x)
}
