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

# TODO: let highlight be a threshold

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
    value <- x[variable,]
    hl <- col[variable,]

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



## panel functions
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

panel.squares <- panel.rectangles

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

panel.blocks <- panel.tiles

### hl is ignored
panel.lines <- function(value, spacing, hl) {
  grid.lines(
    x = seq(length(value)),
    y = value * (1 - spacing),
    default.units = "native"
  )
}


## add cut lines manually to a bertin plot
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
