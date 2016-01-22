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



bertinplot <-
  function(x, order = NULL, highlight = TRUE, options = NULL)
  {
    if(!is.matrix(x))
      stop("Argument 'x' must be a matrix.")
    if(!is.logical(highlight)) 
      stop("Argument 'highlight' must be a logical.")
    
    ## do labels
    if(!is.null(options$xlab)) rownames(x) <- options$xlab
    if(!is.null(options$ylab)) colnames(x) <- options$ylab
    
    ## order
    if(!is.null(order)) x <- permute(x, order)
    
    ## default plot options
    user_options <- options
    options <- list(
      panel.function     = panel.bars, 
      reverse     = FALSE,
      xlab        = NULL,
      ylab        = NULL,
      frame       = FALSE,
      spacing     = 0.2,
      mar         = c(5, 4, 8, 8),
      gp_labels   = gpar(),
      gp_panels   = gpar(),
      shading	    = FALSE,
      newpage     = TRUE,
      pop         = TRUE
    )
    
    ## check and add the plot options
    if(!is.null(user_options) && length(user_options) != 0) {
      o <- pmatch(names(user_options), names(options))
      
      if(any(is.na(o)))
        stop(sprintf(ngettext(length(is.na(o)),
          "Unknown plot option: %s",
          "Unknown plot options: %s"),
          paste(names(user_options)[is.na(o)],
            collapse = " ")))
      
      options[o] <- user_options
    }
    
    ## note: Bertin switched cols and rows for his display!
    if(options$reverse) {
      x <- t(x)
    }
    
    
    ## panel.blocks has no spacing!
    if(identical(options$panel.function, panel.blocks)) {
      options$spacing <- 0
    }
    
    ## scale each variable in x for plotting (between 0 and 1 or -1 and 1)
    ## this can deal with 0s, na, nan, but plots inf as na
    infs <- is.infinite(x)
    infsign <- sign(x[is.infinite(x)])
    
    scalem <-matrix(apply(abs(x), 2, max, na.rm = TRUE),
      byrow= TRUE, ncol=ncol(x), nrow= nrow(x))
    scalem[scalem==0] <- 1	
    
    x <- x/scalem 
    
    if(any(infs)) x[infs] <- infsign
    
    # x <- x/ matrix(apply(abs(x), 2, max, na.rm = TRUE), 
    #    byrow= TRUE, ncol=ncol(x), nrow= nrow(x))
    ## fix division by zero (if all entries in a row are 0)
    #x[is.nan(x)] <- 0
    
    ## highlight
    if(length(highlight) == 1 && highlight) 
      highlight <- x > matrix(colMeans(x, na.rm = TRUE), 
        ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
    else if(length(highlight) == 1 && !highlight)
      highlight <- matrix(FALSE, ncol = ncol(x), nrow = nrow(x))
    else if(any(dim(x) != dim(highlight)))
      stop("Argument 'highlight' has incorrect dimensions.")
    
    ## shading?
    if(options$shading) {
      highlight <- map(x, c(.8,.1))
      highlight[!is.finite(highlight)] <- 1
      highlight <- matrix(grey(highlight), ncol=ncol(x), nrow=nrow(x))
    }
    
    ncol_x  <- ncol(x)
    
    ## clear page
    if(options$newpage) grid.newpage()
    
    ## create outer viewport
    xlim <- c(options$spacing, nrow(x) + 1 - options$spacing)
    pushViewport(plotViewport(margins = options$mar, 
      layout = grid.layout(ncol_x, 1), 
      xscale = xlim, 
      yscale = c(0, ncol_x), 
      default.units = "native", 
      name = "bertin"))
    
    for (variable in 1:ncol_x) { 
      value <- x[, variable]
      hl <- highlight[, variable]
      
      ## handle neg. values
      if(identical(options$panel.function, panel.bars) ||
          identical(options$panel.function, panel.lines)) {
        ylim <- c(min(value,0, na.rm=TRUE),
          max(value,0, na.rm=TRUE) + options$spacing) 
      }else{
        ylim <- c(0,
          max(abs(value),0.1, na.rm=TRUE))
      }
      
      pushViewport(viewport(layout.pos.col = 1, layout.pos.row = variable,
        xscale = xlim, 
        yscale = ylim, 
        default.units = "native", gp = options$gp_panels))
      
      ## call panel function
      options$panel.function(value, options$spacing, hl)
      
      ## do frame
      if(options$frame) grid.rect(x = 1:length(value), 
        width = 1,
        default.units = "native")
      
      upViewport(1)
    }
    
    spacing_corr <- if(options$spacing <= 0) -options$spacing+0.2 else 0
    
    grid.text(rownames(x), x = 1:nrow(x), y = ncol_x + spacing_corr, 
      rot = 90, just = "left",
      default.units= "native", gp = options$gp_labels)
    
    grid.text(rev(colnames(x)), x = 1 + spacing_corr / nrow(x) / 4, 
      y = 0.5:(ncol_x-0.5)/ncol_x,
      just = "left", 
      default.units= "npc", gp = options$gp_labels)
    
    if (options$pop)
      popViewport(1)
    else
      upViewport(1)
  }



## panel functions
panel.bars <- function(value, spacing, hl) {
  grid.rect(x = 1:length(value), y = spacing/2, 
    width = 1 - spacing, 
    height = value*(1 - spacing),
    just = c("centre", "bottom"), default.units = "native", 
    gp = gpar(fill = hl))
}

panel.circles <- function(value, spacing, hl) {
  ## neg. values are dashed
  lty <- as.integer(value<0)+1L
  lty[!is.finite(lty)] <- 0L
  
  value <- abs(value)
  
  grid.circle(x = 1:length(value), y=.5, 
    r = value/2*(1 - spacing),
    default.units = "native", 
    gp = gpar(fill = hl, lty=lty))
}

panel.squares <- function(value, spacing, hl) {
  ## neg. values are dashed
  lty <- as.integer(value<0)+1L
  lty[!is.finite(lty)] <- 0L
  
  grid.rect(x = 1:length(value), 
    width = value*(1 - spacing), 
    height = value*(1 - spacing),
    default.units = "native",
    just = c("centre", "center"),
    gp = gpar(fill = hl, lty=lty))
}

panel.blocks <- function(value, spacing, hl) {
  grid.rect(x = 1:length(value), 
    width = 1, 
    height = unit(1, "npc"),
    default.units = "native",
    just = c("centre", "center"),
    gp = gpar(fill = hl))
}

panel.lines <- function(value, spacing, hl) {
  grid.lines(x = c(1:length(value)), y = value*(1-spacing), 
    default.units = "native")
}


## add cut lines manually to a bertin plot 
bertin_cut_line <- function(x = NULL, y = NULL) {
  if(length(x) <2) x <- rep(x,2)
  if(length(y) <2) y <- rep(y,2)
  
  ## find the bertin Viewport
  if(inherits(try(seekViewport("bertin"), silent=TRUE), "try-error")) {
    stop("bertinplot() needs to be called with options=list(pop=FALSE) first!")
  }
  
  if(is.null(x)) x <- unit(c(0,1), units="npc") 
  else x <- x+.5
  
  if(is.null(y)) y <- unit(c(0,1), units="npc")
  else y <- y
  
  grid.lines(x = x,
    y = y,
    default.units= "native",
    gp=gpar(col="grey", lwd=2))
}
