#######################################################################
# Code to map between ranges for continuous variables
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

## mapping helper

map <- function(x,
  range = c(0, 1),
  from.range = NA) {
  ## deal with infinite values
  infs <- is.infinite(x)
  if (any(infs)) {
    warning(
      "x contains infinite values. +Inf will be mapped to be mapped to largest value + range and -Inf to smallest value - range."
    )
    min_max <- range(x[!infs], na.rm = TRUE)
    pos_inf_val <- min_max[2] + (min_max[2] - min_max[1])
    neg_inf_val <- min_max[1] - (min_max[2] - min_max[1])

    x[infs] <-
      ifelse(sign(x[infs] > 0), pos_inf_val, neg_inf_val)
  }

  ## set from range
  if (any(is.na(from.range)))
    from.range <- range(x, na.rm = TRUE)

  if (length(from.range) != 2L ||
      from.range[1] > from.range[2])
    stop('from.range needs to contain 2 numbers (upper <= lower bound).')
  from.range_width <- from.range[2] - from.range[1]

  if (length(range) != 2L)
    stop('range needs to contain 2 numbers (upper and lower bound).')
  range_width <- range[2] - range[1]

  ## if all values are the same and no from.range is given, then return the average range
  if (from.range_width == 0) {
    x[] <- mean(range)
    return(x)
  }

  ## map to [0,1]
  x <- (x - from.range[1]) / from.range_width

  ## map from [0,1] to [range]
  x <- x * range_width + range[1]

  x
}

map_int <- function(x,
  range = c(1L, 100L),
  from.range = NA) {
  if (length(range) == 1L)
    range <- c(1L, range)
  as.integer(map(x, c(range[1], range[2]), from.range))
}
