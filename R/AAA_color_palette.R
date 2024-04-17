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

#' Different Useful Color Palettes
#'
#' Defines several color palettes for [pimage()], [dissplot()] and
#' [hmap()].
#'
#' The color palettes are created with [colorspace::sequential_hcl()] and
#' [colorspace::diverging_hcl()].
#'
#' The two sequential palettes are: `reds()` and `grays()` (or
#' `greys()`).
#'
#' The two diverging palettes are: `bluered()` and `greenred()`.
#'
#' @name palette
#' @aliases palette, colors
#' @family plots
#'
#' @param n number of different colors produces.
#' @param power used to control how chroma and luminance is increased (1 =
#' linear, 2 = quadratic, etc.)
#' @param bias a positive number. Higher values give more widely spaced colors
#' at the high end.
#' @param ...  further parameters are passed on to [sequential_hcl()]
#' or [diverging_hcl()].
#' @return A vector with `n` colors.
#' @author Michael Hahsler
#' @keywords hplot
#' @examples
#' m <- outer(1:10,1:10)
#' m
#'
#' pimage(m)
#' pimage(m, col = greys(100, power = 2))
#' pimage(m, col = greys(100, bias = 2))
#' pimage(m, col = bluered(100))
#' pimage(m, col = bluered(100, power = .5))
#' pimage(m, col = bluered(100, bias = 2))
#' pimage(m - 25, col = greenred(20, bias = 2))
#'
#' ## choose your own color palettes
#' library(colorspace)
#' hcl_palettes(plot = TRUE)
#'
#' ## blues (with 20 shades)
#' pimage(m,
#'   col = colorspace::sequential_hcl(20, "Blues", rev = TRUE))
#' ## blue to green (aka "Cork")
#' pimage(m,
#'   col = colorspace::diverging_hcl(100, "Cork"))
#' @export
bluered <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(colorspace::diverging_hcl(n, palette = "Blue-Red", power = power, ...),
    bias = bias)(n)

#hclplot(bluered(10))
#plot(1:20, col = bluered(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
greenred <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::diverging_hcl(n, palette = "Red-Green", power = power, ...)
  ), bias = bias)(n)

#hclplot(greenred(10))
#plot(1:20, col = greenred(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
reds <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Reds", power = power, ...)
  ), bias = bias)(n)

#hclplot(reds(10))
#plot(1:20, col = reds(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
blues <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Blues 2", power = power, ...)
  ), bias = bias)(n)

#hclplot(blues(10))
#plot(1:20, col = blues(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
greens <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Greens", power = power, ...)
  ), bias = bias)(n)

#hclplot(greens(10))
#plot(1:20, col = greens(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
greys <- function(n = 100,
  bias = 1,
  power = 1,
  ...)
  grDevices::colorRampPalette(rev(
    colorspace::sequential_hcl(n, palette = "Grays", power = power, ...)
  ), bias = bias)(n)

#hclplot(greys(10))
#plot(1:20, col = greys(20), pch = 19, cex = 4)

#' @rdname palette
#' @export
grays <- greys


.map_color_01 <- function(x, col) {
  x[] <- col[map_int(x, length(col), from.range = c(0, 1))]
  x
}

# translate all data to a color
.map_color <- function(x, col, from.range = NA) {
  x[] <- col[map_int(x, length(col), from.range)]
  x
}

## define default colors
#.sequential_pal <- grays
.sequential_pal <- blues
.diverge_pal <- bluered

## define default ggplot2 colors
.gg_logical_pal <- function()
  ggplot2::scale_fill_manual(values = c("white", "black"), na.value = "white")

.gg_sequential_pal <- function(dist = FALSE, limits = NULL) {
  if (dist)
    ggplot2::scale_fill_gradient(low = scales::muted("blue"),
      high = "white",
      na.value = "white",
      limits = limits)
  else
    ggplot2::scale_fill_gradient(low = "white",
      high = scales::muted("blue"),
      na.value = "white",
      limits = limits)
}

.gg_diverge_pal <- function(limits = NULL)
  ggplot2::scale_fill_gradient2(
    low = scales::muted("blue"),
    mid = "white",
    high = scales::muted("red"),
    na.value = "white",
    midpoint = 0,
    limits = limits
  )
