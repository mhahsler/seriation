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

# brower and kile 1988
concentrate <- function(x){
  
  # step 1: calculate mean column position (mcp) of presences across rows
  mcp <- unlist(apply(
    x, 
    MARGIN = 1, 
    FUN = function(z){ mean(which(z == 1)) }, 
    simplify = FALSE
  ))
  
  # step 2: sort rows by mcp
  x <- x[order(mcp), ]
  
  # step 3: calculate mean row position (mrp) of presences across columns
  mrp <- unlist(apply(
    x, 
    MARGIN = 2, 
    FUN = function(z){ mean(which(z == 1)) }, 
    simplify = FALSE
  ))
  
  # step 4: sort columns by mrp
  x[, order(mrp)]
  
}

seriate_bku <- function(x, control = NULL){

  old <- x
  
  not_identical <- TRUE
  
  while(not_identical){
    
    new <- concentrate(old)
    
    not_identical <- !identical(old, new)
    
    old <- new
    
  }
  
  new

}

set_seriation_method(
  kind = "matrix",
  name = "BK_unconstrained",
  definition = seriate_bku,
  description = "Order using mean row and column position of presences (1s)."
)