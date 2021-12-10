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

# library(colorspace)

.map_color_01 <- function(x, col) {
  x[] <- col[map_int(x, length(col), from.range = c(0, 1))]
  x
}

# translate all data to a color
.map_color <- function(x, col, from.range = NA) {
  x[] <- col[map_int(x, length(col), from.range)]
  x
}

bluered <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(colorspace::diverging_hcl(n, palette = "Blue-Red", power = power, ...),
    bias = bias)(n)

#hclplot(bluered(10))
#plot(1:20, col = bluered(20), pch = 19, cex = 4)

greenred <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::diverging_hcl(n, palette = "Red-Green", power = power, ...)
  ), bias = bias)(n)

#hclplot(greenred(10))
#plot(1:20, col = greenred(20), pch = 19, cex = 4)

reds <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Reds", power = power, ...)
  ), bias = bias)(n)

#hclplot(reds(10))
#plot(1:20, col = reds(20), pch = 19, cex = 4)

blues <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Blues 2", power = power, ...)
  ), bias = bias)(n)

#hclplot(blues(10))
#plot(1:20, col = blues(20), pch = 19, cex = 4)

greens <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Greens", power = power, ...)
  ), bias = bias)(n)

#hclplot(greens(10))
#plot(1:20, col = greens(20), pch = 19, cex = 4)

greys <- grays <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Grays", power = power, ...)
  ), bias = bias)(n)

#hclplot(greys(10))
#plot(1:20, col = greys(20), pch = 19, cex = 4)

## define default colors
#.sequential_pal <- grays
.sequential_pal <- blues
.diverge_pal <- bluered

## define default ggplot2 colors
.gg_logical_pal <- function()
   ggplot2::scale_fill_manual(values = c("white", "black"), na.value = "white")

.gg_sequential_pal <- function(dist = FALSE) {
  if (dist)
    ggplot2::scale_fill_gradient(low = scales::muted("blue"),
      high = "white",
      na.value = "white")
  else
    ggplot2::scale_fill_gradient(low = "white",
      high = scales::muted("blue"),
      na.value = "white")
}

.gg_diverge_pal <- function()
  ggplot2::scale_fill_gradient2(
    low = scales::muted("red"),
    mid = "white",
    high = scales::muted("blue"),
    na.value = "white",
    midpoint = 0
  )
