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
#' type (\code{kind}) (e.g., "dist", "matrix").
#' The result is a vector of character strings with the
#' method names that can be used in function `seriate()`.
#' If \code{kind} is missing, then a list of
#' methods is returned.
#'
#' \code{get_seriation_method()} returns detailed information for a given method in
#' form of an object of class \code{"seriation_method"}.
#' The information includes a description, parameters and the
#' implementing function.
#'
#' With \code{set_seriation_method()} new seriation methods can be added by the
#' user. The implementing function (\code{definition}) needs to have the formal
#' arguments \code{x, control} and, for arrays and matrices \code{margin},
#' where \code{x} is the data object and
#' \code{control} contains a list with additional information for the method
#' passed on from \code{seriate()}, and \code{margin} is a vector specifying
#' what dimensions should be seriated.
#' The implementation has to return a list of
#' objects which can be coerced into \code{ser_permutation_vector} objects
#' (e.g., integer vectors). The elements in the list have to be in
#' corresponding order to the dimensions of \code{x}.
#'
#' @import registry
#' @name registry_for_seriaiton_methods
#' @family seriation
#'
#' @param kind the data type the method works on. For example, \code{"dist"},
#' \code{"matrix"} or \code{"array"}. If missing, then methods for any type are
#' shown.
#' @param name the name for the method used to refer to the method in
#' [seriate()].
#' @param names_only logical; return only the method name. `FALSE` returns
#'    also the method descriptions.
#' @param definition a function containing the method's code.
#' @param description a description of the method. For example, a long name.
#' @param control a list with control arguments and default values.
#' @param randomized logical; does the algorithm use randomization and re-running
#'   the algorithm several times will lead to different results (see: [seriate_rep()]).
#' @param optimizes what criterion does the algorithm try to optimize
#'   (see: [list_criterion_methods()]).
#' @param x an object of class  "seriation_method" to be printed.
#' @param verbose logical; print a message when a new method is registered.
#' @param ... further information that is stored for the method in the
#' registry.
#' @returns
#' - \code{list_seriation_method()} result is a vector of character strings with the
#'   names of the methods. These names are used for methods in `seriate()`.
#' - \code{get_seriation_method()} returns a given method in form of an object of class
#'   \code{"seriation_method"}.
#'
#' @author Michael Hahsler
#' @seealso This registry uses [registry::registry].
#' @keywords misc
#' @examples
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
#' # 1. Create the seriation method: Reverse the row order
#' #    (NA means no seriation is applied to columns)
#' seriation_method_reverse_rows <- function(x, control = NULL, margin = c(1, 2)) {
#'     list(rev(seq(nrow(x))), NA)[margin]
#' }
#'
#' # 2. Register new method
#' set_seriation_method("matrix", "Reverse_rows", seriation_method_reverse_rows,
#'     description = "Reverse identity order", control = list())
#'
#' list_seriation_methods("matrix")
#' get_seriation_method("matrix", "reverse_rows")
#'
#' # 3. Use the new seriation methods
#' seriate(matrix(1:12, ncol = 3), "reverse_rows")
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

registry_seriate$set_field("randomized", type = "logical",
                           is_key = FALSE)

registry_seriate$set_field("optimizes", type = "character",
                           is_key = FALSE)

registry_seriate$set_field("registered_by", type = "character",
                             is_key = FALSE)

#' @rdname registry_for_seriaiton_methods
#' @export
list_seriation_methods <- function(kind, names_only = TRUE) {
  if (missing(kind)) {
    kinds <- unique(sort(as.vector(
      sapply(registry_seriate$get_entries(), "[[", "kind")
    )))

    sapply(
      kinds,
      FUN = function(k)
        list_seriation_methods(k, names_only = names_only)
    )

  } else{
    if (names_only)

      sort(as.vector(sapply(
        registry_seriate$get_entries(kind = kind), "[[", "name"
      )))
    else {
      l <- registry_seriate$get_entries(kind = kind)
      l[order(names(l))]
    }
  }
}

#' @rdname registry_for_seriaiton_methods
#' @export
get_seriation_method <- function(kind, name) {

  ## catch deprecated methods
  if (tolower(name) == "mds_nonmetric") {
    name <- "isoMDS"
    warning("seriation method 'MDS_nonmetric' is now deprecated and will be removed in future releases. Using `isoMDS`")
  }

  if (tolower(name) == "mds_metric") {
    name <- "MDS"
    warning("seriation method 'MDS_metric' is now deprecated and will be removed in future releases. Using `MDS`")
  }


  if (missing(kind)) {
    method <- registry_seriate$get_entry(name = name)
    kind <- NA
  }  else
    method <- registry_seriate$get_entry(kind = kind, name = name)

  if (is.null(method))
    stop(
      "Unknown seriation method ",
      name,
      " for data type ",
      kind,
      ". Maybe the method has not been registered yet. ",
      "Check list_seriation_methods()."
    )

  method
}

#' @rdname registry_for_seriaiton_methods
#' @export
set_seriation_method <- function(kind,
                                 name,
                                 definition,
                                 description = NULL,
                                 control = list(),
                                 randomized = FALSE,
                                 optimizes = NA_character_,
                                 verbose = FALSE,
                                 ...) {
  ## check formals
  if (!identical(names(formals(definition)),
                 c("x", "control")) &&
      !identical(names(formals(definition)),
                 c("x", "control", "margin")))
    stop("Seriation methods must have formals 'x', 'control' and optionally 'margin'.")

  if (sys.nframe() > 1) {
    caller <- deparse(sys.calls()[[sys.nframe()-1]])
    if (is.null(caller) || !startsWith(caller, "register_"))
    caller <- NA_character_
  } else
    caller <- "manual"

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
      name = name,
      kind = kind,
      fun = definition,
      description = description,
      control = control,
      randomized = randomized,
      optimizes = optimizes,
      registered_by = caller
    )
  } else {
    registry_seriate$set_entry(
      name = name,
      kind = kind,
      fun = definition,
      description = description,
      control = control,
      randomized = randomized,
      optimizes = optimizes,
      registered_by = caller
    )
  }

  if (verbose)
    message("Registering new seriation method ",
        sQuote(name),
        " for ",
        sQuote(kind),
        caller
        )
}


#' @rdname registry_for_seriaiton_methods
#' @export
print.seriation_method <- function(x, ...) {
  if (is.na(x$optimizes))
    opt <- "Other"
  else
    opt <- x$optimizes

  if (!is.null(attr(x$optimizes, "description")))
      opt <- paste0(opt, " (", attr(x$optimizes, "description"), ")")

  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    gettextf("optimizes:   %s", opt),
    gettextf("randomized:  %s", x$randomized)
  ))

  if(!is.na(x$registered_by))
    writeLines(gettextf("registered by: %s", x$registered_by))

  writeLines(c(
    strwrap(
      gettextf("description: %s", x$description),
      prefix = "             ",
      initial = ""
    )
  ))

  writeLines("control:")
  .print_control(x$control)

  invisible(x)
}


.print_control <- function(control,
                           label = "default values",
                           help = TRUE,
                           trim_values = 30L) {
  if (length(control) < 1L) {
    writeLines("no parameters")
  } else{
    contr <- lapply(
      control,
      FUN = function(x)
        strtrim(paste(deparse(x), collapse = ""), trim_values)
    )

    contr <- as.data.frame(t(as.data.frame(contr)))
    colnames(contr) <- c(label)

    contr <- cbind(contr, help = "N/A")
    if (!is.null(attr(control, "help")))
      for (i in seq(nrow(contr))) {
        hlp <- attr(control, "help")[[rownames(contr)[i]]]
        if (!is.null(hlp))
        contr[["help"]][i] <- hlp
      }

    print(contr, quote = FALSE)
  }

  cat("\n")
}

.opt <- function(criterion, description = NULL)
  structure(criterion, description = description)
