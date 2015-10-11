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
    

## setup registries

## seriate
registry_seriate <- registry(registry_class="seriation_registry", 
  entry_class="seriation_method")

registry_seriate$set_field("kind", type = "character", 
  is_key = TRUE, index_FUN = match_partial_ignorecase)
registry_seriate$set_field("name", type = "character", 
  is_key = TRUE, index_FUN = match_partial_ignorecase)
registry_seriate$set_field("fun", type = "function", 
  is_key = FALSE)
registry_seriate$set_field("description", type = "character", 
  is_key = FALSE)


print.seriation_method <- function(x, ...) {
  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    gettextf("description: %s", x$description)))
  invisible(x)
}


## criterion
registry_criterion <- registry(registry_class="criterion_registry", 
  entry_class="criterion_method")

registry_criterion$set_field("kind", type = "character", 
  is_key = TRUE, index_FUN = match_partial_ignorecase)
registry_criterion$set_field("name", type = "character", 
  is_key = TRUE, index_FUN = match_partial_ignorecase)
registry_criterion$set_field("fun", type = "function", 
  is_key = FALSE)
registry_criterion$set_field("description", type = "character", 
  is_key = FALSE)
registry_criterion$set_field("merit", type = "logical", 
  is_key = FALSE)


print.criterion_method <-function(x, ...) {
  writeLines(c(
    gettextf("name:        %s", x$name),
    gettextf("kind:        %s", x$kind),
    gettextf("description: %s", x$description),
    gettextf("merit:       %s", x$merit)))
  invisible(x)
}



