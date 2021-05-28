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


bluered <- function(n = 100, bias = 1, power = 1, ...)
  grDevices::colorRampPalette(colorspace::diverging_hcl(n, palette = "Blue-Red", power = power, ...), bias = bias)(n)

greenred <- function(n = 100, bias = 1, power = 1,...)
  grDevices::colorRampPalette(rev(colorspace::diverging_hcl(n, palette = "Red-Green", power = power, ...)), bias = bias)(n)

reds <- function(n = 100, bias = 1, power = 1, ...)
  grDevices::colorRampPalette(rev(colorspace::sequential_hcl(n, palette = "Reds", power = power, ...)), bias = bias)(n)

greys <- grays <- function(n = 100, bias = 1, power = 1, ...)
  grDevices::colorRampPalette(rev(colorspace::sequential_hcl(n, palette = "Grays", power = power, ...)), bias = bias)(n)

## define default colors
.sequential_pal <- grays
.diverge_pal <- bluered
