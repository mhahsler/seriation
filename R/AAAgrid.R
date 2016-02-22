#######################################################################
# Basic Grid helpers
# Copyrigth (C) 2011 Michael Hahsler
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



## grid helpers

.grid_basic_layout <- function(main = "", 
  left=unit(4, "lines"), right = unit(4, "lines"), 
  bottom = unit(4, "lines"),
  gp = gpar()){
  
  pushViewport(viewport(layout = grid.layout(4, 3,
    widths = unit.c(
      left,                             # space
      unit(1, "npc") - left - right,    # plot
      right                             # space
    ),
    heights = unit.c(
      unit(3, "lines"),                  # title
      unit(1, "lines"),                  # space
      unit(1, "npc") - unit(4, "lines") - bottom, # plot
      bottom                            # space
    )
  ), gp = gp))
  
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1, 
    name="main"))

  gp$cex <- 1.3
  gp$fontface <- "bold"
  
  grid.text(main, gp = gp)
  upViewport(1)
  
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3, 
    name="plot"))
  upViewport(2)
}

.grid_basic_layout_with_colorkey <- function(main = "",
  left=unit(4, "lines"), right = unit(4, "lines"), 
  bottom = unit(4, "lines"),
  gp = gpar()){
  
  pushViewport(viewport(layout = grid.layout(4, 3,
    widths = unit.c(
      left,                             # space
      unit(1, "npc") - left - right,    # plot
      right                             # space
    ),
    heights = unit.c(
      unit(3, "lines"),                  # title
      unit(1, "lines"),                  # space
      unit(1, "npc") - unit(4, "lines") - bottom, # plot
      bottom                            # space
    )
  ), gp = gp))
  
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1, 
    name="main"))
  
  gp$cex <- 1.3
  gp$fontface <- "bold"
  grid.text(main, gp = gp)
  upViewport(1)
  
  pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 3))
  
  pushViewport(viewport(layout = grid.layout(1, 3,
    widths = unit.c(
      unit(1, "npc") - unit(8, "lines"),# plot
      unit(1, "lines"),                  # space
      unit(1, "lines")                  # colorkey
    ),
    heights = unit.c(
      unit(1, "npc") # plot
    )
  )))
  
  pushViewport(viewport(layout.pos.col = 1, layout.pos.row = 1, 
    name="plot"))
  
  upViewport(1)
  
  pushViewport(viewport(layout.pos.col = 3, layout.pos.row = 1, 
    name="colorkey"))
  
  upViewport(1)
  upViewport(2)
}



### new version below uses grid.raster
.grid_image_old <- function(x, y, z, zlim, col = grey.colors(12, 1, 0), 
  name = "image", gp = gpar()) {
  
  if(is.matrix(x)){ 
    z <- x
    x <- 1:ncol(z)
    y <- 1:nrow(z)
  }
  
  if(missing(zlim)) zlim <- range(z, na.rm = TRUE)
  else {# fix data for limits
    z[z < zlim[1]] <- NA
    z[z > zlim[2]] <- NA
  }
  
  offset <- if(zlim[1] < 0) -zlim[1] else 0
  range <- diff(zlim) 
  
  div <- 1/length(col)
  
  ## create a viewport
  vp <- viewport(
    xscale = c(0,(length(x)+1)), yscale = c((length(y)+1),0),
    default.units = "native", name = name)
  pushViewport(vp)
  
  ## make sure we have a color for the maximal value (see floor +1)
  col[length(col)+1] <- col[length(col)]
  
  ## the highest value is lightest color!
  xs <- sapply(x, "rep.int", times = length(y))
  grid.rect(x = xs, y = y, 1, 1, 
    gp = gpar(fill = col[floor((z + offset)/range/div)+1], col=0), 
    default.units = "native")
  
  ## make border
  gp_border       <- gp
  gp_border$fill  <- "transparent"
  grid.rect(x = (length(x)+1)/2, y = (length(y)+1)/2, 
    width = length(x), height = length(y), 
    default.units = "native", gp = gp_border)
  
  upViewport(1)
}

.grid_image <- function(x, zlim, col = grey.colors(12, 1, 0), 
  prop = FALSE, name = "image", gp = gpar()) {
  
  if(missing(zlim)) zlim <- range(x, na.rm = TRUE)
  else {# fix data for limits
    x[x < zlim[1]] <- NA
    x[x > zlim[2]] <- NA
  }
  
  offset <- if(zlim[1] < 0) -zlim[1] else 0
  #range <- diff(zlim) ### not al plots have 0 as the minimum! 
  range <- zlim[2]+offset 
  
  div <- 1/length(col)
  
  ## create a viewport
  if(!prop) {
    vp <- viewport(
      #xscale = c(0,ncol(x)), yscale = c(nrow(x),0),
      xscale = c(0.5,ncol(x)+.5), yscale = c(nrow(x)+.5,0.5),
      default.units = "native", name = name)
    pushViewport(vp)
  }else{
    ## ratio
    if (nrow(x) > ncol(x)) { w <- ncol(x)/nrow(x); h <- 1 } 
    else if (nrow(x) < ncol(x)) { h <- nrow(x)/ncol(x); w <- 1 }
    else { w <- 1; h <- 1 }
    vp <- viewport(
      #xscale = c(0,ncol(x)), yscale = c(nrow(x),0),
      xscale = c(0.5,ncol(x)+.5), yscale = c(nrow(x)+.5,0.5),
      width = unit(w, "snpc"), height = unit(h, "snpc"),
      default.units = "native", name = name)
    pushViewport(vp)
  }
  
  ## make sure we have a color for the maximal value (see floor +1)
  col[length(col)+1] <- col[length(col)]
  
  ## the highest value is lightest color!
  x[] <- col[floor((x + offset)/range/div)+1]
  grid.raster(x, interpolate=FALSE, default.units = "npc", width=1, height=1)
  
  ## make border
  gp_border       <- gp
  gp_border$fill  <- "transparent"
  grid.rect(gp = gp_border)
  
  upViewport(1)
}




.grid_barplot_horiz <- function(height, name = "barplot", xlab="", 
  gp = gpar(), gp_bars = gpar(fill="lightgrey")) {
  
  n <-  length(height)
  
  ## these plots always start at x = 0 or below!
  lim <- c(min(c(height, 0)), max(height))
  
  ## create a viewport
  vp <- viewport(
    xscale = lim , yscale = c(n,0), default.units = "native", name = name, 
    gp =gp)
  pushViewport(vp)
  
  grid.rect(x = 0, y = (1:n)-.5, width = height, height = 1,
    just = c("left", "center"), default.units = "native", gp =gp_bars)
  
  ## hopefully there is space outside for axes
  grid.xaxis()
  grid.text(xlab, y = unit(-3, "lines"))
  
  upViewport(1)
}

# new colorkey uses grid.raster
.grid_colorkey_old <- function(range, col, threshold = NULL, 
  name = "colorkey", gp = gpar()) {
  
  vp <- viewport(
    xscale = range, yscale = c(0,1), 
    default.units = "native", name = name)
  pushViewport(vp)
  
  n <- length(col) 
  width <- diff(range)/n
  xs <- seq(range[1] + width/2, range[2] - width/2, length.out = n)
  
  ## do not display the part above the threshold 
  col[xs > threshold] <- NA
  
  ## col
  gp_col      <- gp
  gp_col$col  <- 0
  gp_col$fill <- col
  grid.rect(x = xs, y = 0, width = width, height = 1,
    just = c("centre", "bottom"), default.units = "native",
    gp = gp_col)
  
  
  ## box
  gp_border       <- gp
  gp_border$fill  <- "transparent"
  grid.rect(x = 0, y = 0, width = 1, height = 1,
    just = c("left", "bottom"), default.units = "npc",
    gp = gp_border)
  
  grid.xaxis(gp = gp)
  
  upViewport(1)
}

.grid_colorkey <- function(range, col, threshold = NULL, lab = "",
  name = "colorkey", horizontal=TRUE, gp = gpar()) {
  
  if(horizontal)
    vp <- viewport(
      xscale = range, yscale = c(0,1), 
      default.units = "native", name = name)
  else	
    vp <- viewport(
      xscale = c(0,1), yscale = range, 
      default.units = "native", name = name)
  
  pushViewport(vp)
  
  n <- length(col) 
  #width <- diff(range)/n
  #xs <- seq(range[1] + width/2, range[2] - width/2, length.out = n)
  xs <- seq(range[1], range[2], length.out = n)
  
  ## do not display the part above the threshold 
  col[xs > threshold] <- NA
  
  ## col
  if(horizontal) grid.raster(t(col), width=1, height=1, interpolate=FALSE)
  else grid.raster(rev(col), width=1, height=1, interpolate=FALSE)
  
  #gp_col      <- gp
  #gp_col$col  <- 0
  #gp_col$fill <- col
  #grid.rect(x = xs, y = 0, width = width, height = 1,
  #    just = c("centre", "bottom"), default.units = "native",
  #    gp = gp_col)
  
  
  ## box
  gp_border       <- gp
  gp_border$fill  <- "transparent"
  grid.rect(gp = gp_border)
  
  if(horizontal) grid.xaxis(gp = gp)
  else grid.yaxis(main = FALSE, gp = gp)
  
  if(horizontal) grid.text(lab, y= unit(-2.5, "lines"))  
  else grid.text(lab, x= unit(3.5, "lines"), rot = 90)
  
  upViewport(1)
}
