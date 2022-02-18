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

#' Plot a Bertin Matrix
#'
#' Plot a data matrix of cases and variables. Each value is represented by a
#' symbol. Large values are highlighted. Note that Bertin arranges the cases
#' horizontally and the variables as rows. The matrix can be rearranged using
#' seriation techniques to make structure in the data visible (see Falguerolles
#' et al 1997).
#'
#' The plot is organized as a matrix of symbols. The symbols are drawn by a
#' panel function, where all symbols of a row are drawn by one call of the
#' function (using vectorization). The interface for the panel function is
#' `panel.myfunction(value, spacing, hl)`. `value` is the vector of
#' values for a row scaled between 0 and 1, `spacing` contains the
#' relative space between symbols and `hl` is a logical vector indicating
#' which symbol should be highlighted.
#'
#' Cut lines can be added to an existing Bertin plot using
#' `bertin_cut_line(x = NULL, y = NULL)`. `x`/`y` is can be a
#' number indicating where to draw the cut line between two columns/rows. If
#' both `x` and `y` is specified then one can select a row/column and
#' the other can select a range to draw a line which does only span a part of
#' the row/column. It is important to call `bertinplot()` with the option
#' `pop = FALSE`.
#'
#' `ggbertinplot()` calls [ggpimage()] and all additional parameters are
#' passed on.
#'
#' @family plots
#' @param x a data matrix. Note that following Bertin, columns are variables
#' and rows are cases. This behavior can be reversed using `reverse = TRUE`
#' in `options`.
#' @param order an object of class `ser_permutation` to rearrange `x`
#' before plotting.  If `NULL`, no rearrangement is performed.
#' @param panel.function a function to produce the symbols. Currently available
#' functions are `panel.bars` (default), `panel.circles`,
#' `panel.rectangles`, `panel.tiles` and `panel.lines`. For
#' circles and squares neg. values are represented by a dashed border. For
#' blocks all blocks are the same size (can be used with `shading = TRUE`).
#' @param geom visualization type. Available ggplot2 geometries are: `"tile"`,
#' `"rectangle"`, `"circle"`, `"line"`, `"bar"`, `"none"`.
#' @param highlight a logical scalar indicating whether to use highlighting.
#' If `TRUE`, all variables with values greater than the variable-wise
#' mean are highlighted. To control highlighting, also a logical matrix or a
#' matrix with colors with the same dimensions as `x` can be supplied.
#' @param row_labels,col_labels a logical indicating if row and column labels
#' in `x` should be displayed.  If `NULL` then labels are displayed
#' if the `x` contains the appropriate dimname and the number of labels is
#' 25 or less. A character vector of the appropriate length with labels can
#' also be supplied.
#' @param flip_axes logical indicating whether to swap cases and variables in
#' the plot. The default (`TRUE`) is to plot cases as columns and
#' variables as rows.
#' @param prop logical; change the aspect ratio so cells in the image have a
#' equal width and height.
#' @param col,y and x in `bertin_cut_line()` are for adding a line to a `bertinplot()` (not ggplot2-based).
#' @param value,spacing,hl are used internally for the panel functions.
#' @param ...
#'   `ggbertinplot()`: further parameters are passed on to [ggpimage()].
#'
#'   `bertinplot()`: further parameters can include:
#'  - `xlab, ylab` labels (default: use labels from `x`).
#'  - `spacing` relative space between symbols (default: 0.2).
#'  - `shading` use gray shades to encode value instead of
#'    highlighting (default: `FALSE`).
#'  - `shading.function` a function that accepts a single argument in range \eqn{[.1, .8]}
#'    and returns a valid corresponding color (e.g., using [rgb()]).
#'  - `frame` plot a grid to separate symbols (default: `FALSE`).
#'  - `mar` margins (see [par()]).
#'  - `gp_labels` `gpar` object for labels (see [gpar()])
#'  - `gp_panels` `gpar` object for panels (see [gpar()]).
#'  - `newpage` a logical indicating whether to start
#'     the plot on a new page (see [grid.newpage()]).
#'  -  `pop` a logical indicating whether to pop the created viewports
#'     (see [pop.viewport()])?
#'
#' @returns Nothing.
#'
#' @author Michael Hahsler
#' @references de Falguerolles, A., Friedrich, F., Sawitzki, G. (1997): A
#' Tribute to J. Bertin's Graphical Data Analysis. In: Proceedings of the
#' SoftStat '97 (Advances in Statistical Software 6), 11--20.
#' @keywords hplot cluster
#' @examples
#' data("Irish")
#' scale_by_rank <- function(x) apply(x, 2, rank)
#' x <- scale_by_rank(Irish[,-6])
#'
#' # Use the the sum of absolute rank differences
#' order <- c(
#'   seriate(dist(x, "minkowski", p = 1)),
#'   seriate(dist(t(x), "minkowski", p = 1))
#' )
#'
#' # Plot
#' bertinplot(x, order)
#'
#' # Some alternative displays
#' bertinplot(x, order, panel = panel.tiles, shading_col = bluered(100), highlight = FALSE)
#' bertinplot(x, order, panel = panel.circles, spacing = -.2)
#' bertinplot(x, order, panel = panel.rectangles)
#' bertinplot(x, order, panel = panel.lines)
#'
#' # Plot with cut lines (we manually set the order here)
#' order <- ser_permutation(c(21, 16, 19, 18, 14, 12, 20, 15,
#'     17, 26, 13, 41,  7, 11, 5, 23, 28, 34, 31, 1, 38, 40,
#'     3, 39,  4, 27, 24,  8, 37, 36, 25, 30, 33, 35,  2,
#'     22, 32, 29, 10,  6,  9),
#'     c(4, 2, 1, 6, 8, 7, 5, 3))
#'
#' bertinplot(x, order, pop=FALSE)
#' bertin_cut_line(, 4) ## horizontal line between rows 4 and 5
#' bertin_cut_line(, 7) ## separate "Right to Life" from the rest
#' bertin_cut_line(14, c(0, 4)) ## separate a block of large values (vertically)
#'
#' # ggplot2-based plots
#' if (require("ggplot2")) {
#'   library(ggplot2)
#'
#'   # Default plot uses bars and highlighting values larger than the mean
#'   ggbertinplot(x, order)
#'
#'   # highlight values in the 4th quartile
#'   ggbertinplot(x, order, highlight = quantile(x, probs = .75))
#'
#'   # Use different geoms. "none" lets the user specify their own geom.
#'   # Variables set are row, col and x (for the value).
#'
#'   ggbertinplot(x, order, geom = "tile", prop = TRUE)
#'   ggbertinplot(x, order, geom = "rectangle")
#'   ggbertinplot(x, order, geom = "rectangle", prop = TRUE)
#'   ggbertinplot(x, order, geom = "circle")
#'   ggbertinplot(x, order, geom = "line")
#'
#'   # Tiles with diverging color scale
#'   ggbertinplot(x, order, geom = "tile", prop = TRUE) +
#'     scale_fill_gradient2(midpoint = mean(x))
#'
#'   # Custom geom (geom = "none"). Defined variables are row, col, and x for the value
#'   ggbertinplot(x, order, geom = "none", prop = FALSE) +
#'     geom_point(aes(x = col, y = row, size = x, color = x > 30), pch = 15) +
#'     scale_size(range = c(1, 10))
#'
#'   # Use a ggplot2 theme with theme_set()
#'   old_theme <- theme_set(theme_minimal() +
#'       theme(panel.grid = element_blank())
#'     )
#'   ggbertinplot(x, order, geom = "bar")
#'   theme_set(old_theme)
#' }
#' @export
bertinplot <- function(x,
  order = NULL,
  panel.function = panel.bars,
  highlight = TRUE,
  row_labels = TRUE,
  col_labels = TRUE,
  flip_axes = TRUE,
  ...) {
  if (!is.matrix(x))
    stop("Argument 'x' must be a matrix.")

  # add ... to options
  options <- list(...)
  options$panel.function <- panel.function

  options <- .get_parameters(
    options,
    list(
      panel.function = panel.bars,
      flip_axes   = TRUE,
      frame       = FALSE,
      spacing     = 0.2,
      margins     = c(5, 4, 8, 8),
      gp_labels   = gpar(),
      gp_panels   = gpar(),
      shading	    = NULL,
      shading_col = .sequential_pal(100),
      newpage     = TRUE,
      pop         = TRUE
    )
  )

  ## panel.blocks has no spacing!
  if (identical(options$panel.function, panel.blocks))
    options$spacing <- 0

  if (is.null(options$shading))
    if (identical(options$panel.function, panel.blocks)) {
      options$shading <- TRUE
    } else {
      options$shading <- FALSE
    }

  ## order
  if (!is.null(order))
    x <- permute(x, order)

  ## note: Bertin switched cols and rows for his display!
  # change x and y?
  if (flip_axes) {
    x <- t(x)
    tmp <- row_labels
    row_labels <- col_labels
    col_labels <- tmp
  }

  ## highlight
  if (is.logical(highlight) && highlight)
    highlight <- mean(x, na.rm = TRUE)

  ## clear page
  if (options$newpage)
    grid.newpage()

  ## create outer viewport
  xlim <- c(options$spacing, ncol(x) + 1 - options$spacing)
  pushViewport(
    plotViewport(
      margins = options$mar,
      layout = grid.layout(nrow(x), 1),
      xscale = xlim,
      yscale = c(0, nrow(x)),
      default.units = "native",
      name = "bertin"
    )
  )

  # shading and highlighting
  if (options$shading)
    col <- .map_color(x, options$shading_col)
  else
    col <- matrix(1, nrow = nrow(x), ncol = ncol(x))

  if (highlight)
    col[x < highlight] <- NA

  # map to [0, 1]
  x <- map(x)

  for (variable in seq(nrow(x))) {
    value <- x[variable, ]
    hl <- col[variable, ]

    ## handle neg. values
    if (identical(options$panel.function, panel.bars) ||
        identical(options$panel.function, panel.lines)) {
      ylim <- c(min(value, 0, na.rm = TRUE),
        max(value, 0, na.rm = TRUE) + options$spacing)
    } else{
      ylim <- c(0,
        max(abs(value), 0.1, na.rm = TRUE))
    }

    pushViewport(
      viewport(
        layout.pos.col = 1,
        layout.pos.row = variable,
        xscale = xlim,
        yscale = ylim,
        default.units = "native",
        gp = options$gp_panels
      )
    )

    ## call panel function
    options$panel.function(value, options$spacing, hl)

    ## do frame
    if (options$frame)
      grid.rect(
        x = seq(length(value)),
        width = 1,
        default.units = "native",
        gp = gpar(fill = NA)
      )

    upViewport(1)
  }

  spacing_corr <-
    if (options$spacing <= 0)
      - options$spacing + 0.2
  else
    0

  grid.text(
    colnames(x),
    x = seq(ncol(x)),
    y = nrow(x) + spacing_corr,
    rot = 90,
    just = "left",
    default.units = "native",
    gp = options$gp_labels
  )

  grid.text(
    rev(rownames(x)),
    x = 1 + spacing_corr / ncol(x) / 4,
    y = 0.5:(nrow(x) - 0.5) / nrow(x),
    just = "left",
    default.units = "npc",
    gp = options$gp_labels
  )

  if (options$pop)
    popViewport(1)
  else
    upViewport(1)
}

#' @rdname bertinplot
#' @export
panel.bars <- function(value, spacing, hl) {
  grid.rect(
    x = seq(length(value)),
    y = spacing / 2,
    width = 1 - spacing,
    height = value * (1 - spacing),
    just = c("centre", "bottom"),
    default.units = "native",
    gp = gpar(fill = hl)
  )
}

#' @rdname bertinplot
#' @export
panel.circles <- function(value, spacing, hl) {
  ## neg. values are dashed
  lty <- as.integer(value < 0) + 1L
  lty[!is.finite(lty)] <- 0L

  value <- abs(value)

  value[value == 0] <- NA ### hide empty squares

  grid.circle(
    x = seq(length(value)),
    y = unit(.5, "npc"),
    r = value / 2 * (1 - spacing),
    default.units = "native",
    gp = gpar(fill = hl, lty = lty)
  )
}

#' @rdname bertinplot
#' @export
panel.rectangles <- function(value, spacing, hl) {
  ## neg. values are dashed
  lty <- as.integer(value < 0) + 1L
  lty[!is.finite(lty)] <- 0L

  value[value == 0] <- NA ### hide emply squares

  grid.rect(
    x = seq(length(value)),
    width = value * (1 - spacing),
    height = value * (1 - spacing),
    default.units = "native",
    just = c("centre", "center"),
    gp = gpar(fill = hl, lty = lty)
  )
}

#' @rdname bertinplot
#' @export
panel.squares <- panel.rectangles

#' @rdname bertinplot
#' @export
panel.tiles <- function(value, spacing, hl) {
  grid.rect(
    x = seq(length(value)),
    width = 1,
    height = unit(1, "npc"),
    default.units = "native",
    just = c("centre", "center"),
    gp = gpar(fill = hl)
  )
}

#' @rdname bertinplot
#' @export
panel.blocks <- panel.tiles

### hl is ignored
#' @rdname bertinplot
#' @export
panel.lines <- function(value, spacing, hl) {
  grid.lines(
    x = seq(length(value)),
    y = value * (1 - spacing),
    default.units = "native"
  )
}


## add cut lines manually to a bertin plot
#' @rdname bertinplot
#' @export
bertin_cut_line <- function(x = NULL,
  y = NULL,
  col = "red") {
  if (length(x) < 2)
    x <- rep(x, 2)
  if (length(y) < 2)
    y <- rep(y, 2)

  ## find the bertin Viewport
  if (inherits(try(seekViewport("bertin"), silent = TRUE)
    , "try-error")) {
    stop("bertinplot() needs to be called with options = list(pop = FALSE) first!")
  }

  if (is.null(x))
    x <- unit(c(0, 1), units = "npc")
  else
    x <- x + .5

  if (is.null(y))
    y <- unit(c(0, 1), units = "npc")
  else
    y <- y

  grid.lines(
    x = x,
    y = y,
    default.units = "native",
    gp = gpar(col = col, lwd = 2)
  )
}
