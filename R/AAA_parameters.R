#######################################################################
# Code to check parameter/control objects
# Copyright (C) 2011 Michael Hahsler
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


## helper to parse parameter lists with defaults
.nodots <- function(...) {
  l <- list(...)
  if (length(l) > 0L)
    warning("Unknown arguments: ",
      paste(names(l), "=", l, collapse = ", "))
}

.get_parameters <- function(parameter, defaults) {
  defaults <- as.list(defaults)
  parameter <- as.list(parameter)

  ## add verbose
  if (is.null(defaults$verbose))
    defaults$verbose <- FALSE

  o <- integer()
  if (length(parameter) != 0) {
    o <- pmatch(names(parameter), names(defaults))

    ## unknown parameter
    if (any(is.na(o))) {
      warning(sprintf(
        "Unknown control parameter(s): %s",
        paste(names(parameter)[is.na(o)],
              collapse = ", ")
      ),
        call. = FALSE,
        immediate. = TRUE)
    }

    ### defaults are now the actual parameters
    defaults[o[!is.na(o)]] <- parameter[!is.na(o)]
  }

  if (defaults$verbose || any(is.na(o))) {
    cat("control:\n")
    .print_control(defaults, "used values")
  }

  defaults
}
