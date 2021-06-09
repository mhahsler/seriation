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
    zlim = NULL,
    key = TRUE,
    keylab = "",
    symkey = TRUE,
    upper_tri = TRUE,
    lower_tri = TRUE,
    row_labels = NULL,
    col_labels = NULL,
    prop = TRUE,
    flip_axes = FALSE,
    reverse_columns = FALSE,
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
    zlim = NULL,
    key = TRUE,
    keylab = "",
    symkey = TRUE,
    upper_tri = TRUE,
    lower_tri = TRUE,
    row_labels = NULL,
    col_labels = NULL,
    prop = TRUE,
    flip_axes = FALSE,
    reverse_columns = FALSE,
    ...,
    newpage = TRUE,
    pop = TRUE,
    gp = NULL) {
    x <- as.matrix(x)

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
    if (any(!upper_tri ||
        !lower_tri) &&
        nrow(x) != ncol(x))
      stop("Upper or lower triangle can only be suppressed for square matrices!")
    if (!upper_tri)
      x[upper.tri(x)] <- NA
    if (!lower_tri)
      x[lower.tri(x)] <- NA

    # change x and y
    if (flip_axes) {
      x <- t(x)
      tmp <- row_labels
      row_labels <- col_labels
      col_labels <- tmp
    }

    # reverse order of columns
    if (reverse_columns)
      x <- x[, seq(ncol(x), 1)]

    # deal with row/col labels
    if (!is.null(row_labels) && !is.logical(row_labels)) {
      if (length(row_labels) != nrow(x))
        stop("Length of row_labels does not match the number of rows of x.")
      rownames(x) <- row_labels
      row_labels <- TRUE
    }

    if (!is.null(col_labels) && !is.logical(col_labels)) {
      if (length(col_labels) != ncol(x))
        stop("Length of col_labels does not match the number of columns of x.")
      colnames(x) <- col_labels
      col_labels <- TRUE
    }

    if (is.null(row_labels))
      if (!is.null(rownames(x)) &&
          nrow(x) < 25) {
        row_labels <- TRUE
      } else{
        row_labels <- FALSE
      }
    if (is.null(col_labels))
      if (!is.null(colnames(x)) &&
          ncol(x) < 25) {
        col_labels <- TRUE
      } else{
        col_labels <- FALSE
      }

    if (is.null(rownames(x)))
      rownames(x) <- seq(nrow(x))
    if (is.null(colnames(x)))
      colnames(x) <- seq(ncol(x))

    # create layout for plot
    bottom_mar <- if (col_labels)
      max(stringWidth(colnames(x))) + unit(3, "lines")
    else
      unit(1, "lines")

    left_mar <- if (row_labels)
      max(stringWidth(rownames(x))) + unit(3, "lines")
    else
      unit(1, "lines")

    if (newpage)
      grid.newpage()

    if (key) {
      .grid_basic_layout_with_colorkey(
        main = main,
        left = left_mar,
        right = unit(0, "lines"),
        bottom = bottom_mar,
        gp = gp
      )
      downViewport("colorkey")
      .grid_colorkey(zlim,
        col = col,
        horizontal = FALSE,
        lab = keylab)
      upViewport(1)

    } else
      .grid_basic_layout(
        main = main,
        left = left_mar,
        right = unit(0, "lines"),
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
    if (col_labels)
      grid.text(
        colnames(x),
        y = unit(-1, "lines"),
        x = unit(1:ncol(x), "native"),
        rot = 90,
        just = "right"
      ) #, gp=gp)
    #grid.xaxis(at=1:ncol(x),
    #	    label=colnames(x))
    if (row_labels)
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
    zlim = NULL,
    key = TRUE,
    keylab = "",
    symkey = TRUE,
    upper_tri = FALSE,
    lower_tri = TRUE,
    row_labels = NULL,
    col_labels = NULL,
    prop = TRUE,
    flip_axes = FALSE,
    reverse_columns = FALSE,
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

    if (flip_axes) warning("flip_axes has no effect for distance matrices.")

    pimage.matrix(
      x,
      order = NULL,  # already reordered
      main = main,
      xlab = xlab,
      ylab = ylab,
      col = col,
      zlim = zlim,
      key = key,
      keylab = keylab,
      symkey = symkey,
      upper_tri = upper_tri,
      lower_tri = lower_tri,
      row_labels = row_labels,
      col_labels = col_labels,
      prop = prop,
      flip_axes = FALSE,
      reverse_columns = reverse_columns,
      ...,
      newpage = newpage,
      pop = pop,
      gp = gp
    )
  }
