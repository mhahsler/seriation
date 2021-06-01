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

## image method that makes a proper image plot of a matrix.
## the rows and columns are swapped and the order of the
## columns (original rows) is reversed.

pimage <-
  function(x,
    order = NULL,
    col = NULL,
    main = "",
    xlab = "",
    ylab = "",
    axes = NULL,
    zlim = NULL,
    key = TRUE,
    key.lab = "",
    symkey = TRUE,
    upper.tri = TRUE,
    lower.tri = TRUE,
    labRow = NULL,
    labCol = NULL,
    prop = NULL,
    ...,
    newpage = TRUE,
    pop = TRUE,
    gp = NULL)
UseMethod("pimage")

### Note for matrix large values are dark, for dist large values are light!
pimage.matrix <-
  function(x,
    order = NULL,
    col = NULL,
    main = "",
    xlab = "",
    ylab = "",
    axes = NULL,
    zlim = NULL,
    key = TRUE,
    key.lab = "",
    symkey = TRUE,
    upper.tri = TRUE,
    lower.tri = TRUE,
    labRow= NULL,
    labCol = NULL,
    prop = NULL,
    ...,
    newpage = TRUE,
    pop = TRUE,
    gp = NULL) {
    x <- as.matrix(x)

    # deprecated
    if(!is.null(axes)) warning("Parameter axes is deprecated and ignored. Use labRow and labCol instead.")

    # check data
    if (all(is.na(x)))
      stop("all data missing in x.")
    if (any(is.infinite(x)))
      stop("x contains infinite entries.")

    # set default values
    # no key for logical data!
    if (is.logical(x))
      key <- FALSE

    if (is.null(col)) {
      if (is.logical(x))
        col <- c("white", "black")
      else if (any(x < 0, na.rm = TRUE))  {
        col <- .diverge_pal(100)
        if (is.null(zlim) && symkey)
          zlim <- max(abs(range(x, na.rm = TRUE))) * c(-1, 1)
      }
      else
        col <- .sequential_pal(100)
    }

    if (is.null(prop))
      prop <- FALSE

    if (is.null(gp))
      gp <- gpar()

    if (is.null(zlim))
      zlim <- range(x, na.rm = TRUE)

    # reorder
    if (!is.null(order))
      x <- permute(x, order)

    # mask triangles
    if (any(!upper.tri ||
        !lower.tri) &&
        nrow(x) != ncol(x))
      stop("Upper or lower triangle can only be suppressed for square matrices!")
    if (!upper.tri)
      x[upper.tri(x)] <- NA
    if (!lower.tri)
      x[lower.tri(x)] <- NA

    # deal with row/col labels
    if (!is.null(labRow) && !is.logical(labRow)) {
      if (length(labRow) != nrow(x))
        stop("Length of labRow does not match the number of rows of x.")
      rownames(x) <- labRow
      labRow <- TRUE
    }

    if (!is.null(labCol) && !is.logical(labCol)) {
      if (length(labCol) != ncol(x))
        stop("Length of labCol does not match the number of columns of x.")
      colnames(x) <- labCol
      labCol <- TRUE
    }

    if (is.null(labRow))
      if (!is.null(rownames(x)) &&
          nrow(x) < 25) {
        labRow <- TRUE
      } else{
        labRow <- FALSE
      }
    if (is.null(labCol))
      if (!is.null(colnames(x)) &&
          ncol(x) < 25) {
        labCol <- TRUE
      } else{
        labCol <- FALSE
      }

    if (is.null(rownames(x)))
      rownames(x) <- seq(nrow(x))
    if (is.null(colnames(x)))
      colnames(x) <- seq(ncol(x))

    # create layout for plot
    bottom_mar <- if (labCol)
      max(stringWidth(colnames(x))) + unit(3, "lines")
    else
      unit(4, "lines")
    left_mar <- if (labRow)
      max(stringWidth(rownames(x))) + unit(3, "lines")
    else
      unit(4, "lines")

    if (newpage)
      grid.newpage()

    if (key) {
      .grid_basic_layout_with_colorkey(
        main = main,
        left = left_mar,
        bottom = bottom_mar,
        gp = gp
      )
      downViewport("colorkey")
      .grid_colorkey(zlim,
        col = col,
        horizontal = FALSE,
        lab = key.lab) #, gp=gp)
      upViewport(1)

    } else
      .grid_basic_layout(
        main = main,
        left = left_mar,
        bottom = bottom_mar,
        gp = gp
      )

    downViewport("plot")
    .grid_image(
      x,
      col = col,
      zlim = zlim,
      prop = prop
    ) #, gp=gp)

    ## axes and labs
    downViewport("image")
    if (labCol)
      grid.text(
        colnames(x),
        y = unit(-1, "lines"),
        x = unit(1:ncol(x), "native"),
        rot = 90,
        just = "right"
      ) #, gp=gp)
    #grid.xaxis(at=1:ncol(x),
    #	    label=colnames(x))
    if (labRow)
      grid.text(
        rownames(x),
        x = unit(-1, "lines"),
        y = unit(1:nrow(x), "native"),
        just = "right"
      ) #, gp=gp)
    #grid.yaxis(at=1:nrow(x),
    #    label=rownames(x))


    if (xlab != "")
      grid.text(xlab, y = -1 * bottom_mar + unit(1, "lines"))
    #, gp=gp)
    if (ylab != "")
      grid.text(ylab,
        x = ,-1 * left_mar + unit(1, "lines"),
        rot = 90) #, gp=gp)

    if (pop)
      popViewport(3)
    else
      upViewport(3)
  }

pimage.default <- pimage.matrix

## small values are dark
pimage.dist <-
  function(x,
    order = NULL,
    col = NULL,
    main = "",
    xlab = "",
    ylab = "",
    axes = NULL,
    zlim = NULL,
    key = TRUE,
    key.lab = "",
    symkey = TRUE,
    upper.tri = TRUE,
    lower.tri = TRUE,
    labRow = NULL,
    labCol = NULL,
    prop = NULL,
    ...,
    newpage = TRUE,
    pop = TRUE,
    gp = NULL) {
    if (is.null(col))
      col <- rev(.sequential_pal(100))
    else
      col <- rev(col)
    if (is.null(prop))
      prop <- TRUE
    if (!is.null(order))
      x <- permute(x, order)

    pimage.matrix(
      x,
      order = NULL,
      main = main,
      xlab = xlab,
      ylab = ylab,
      col = col,
      axes = axes,
      zlim = zlim,
      key = key,
      key.lab = key.lab,
      symkey = symkey,
      upper.tri = upper.tri,
      lower.tri = lower.tri,
      labRow = labRow,
      labCol = labCol,
      prop = prop,
      ...,
      newpage = newpage,
      pop = pop,
      gp = gp
    )
  }
