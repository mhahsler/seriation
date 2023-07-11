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



#' Permutation Image Plot
#'
#' Provides methods for matrix shading, i.e., displaying a color image for
#' matrix (including correlation matrices) and \code{dist} objects given an
#' optional permutation.  The plot arranges colored rectangles to represent the
#' matrix value. Columns and rows appear in the order in the matrix.  This
#' visualization is also know asi a heatmap.  Implementations based on the
#' \pkg{grid} graphics engine and based n \pkg{ggplot2} are provided.
#'
#' Plots a matrix in its original row and column orientation. This means, in a
#' plot the columns become the x-coordinates and the rows the y-coordinates (in
#' reverse order).
#'
#' If \code{x} is of class \code{dist} it is converted to full-storage
#' representation before plotting.
#'
#' \bold{Grid-based plot:} The viewports used for plotting are called:
#' \code{"plot"}, \code{"image"} and \code{"colorkey"}.  \emph{Note:} Since
#' \code{pimage} uses \pkg{grid}, it should not be mixed with base R primitive
#' plotting functions, but the appropriate functions in
#' \code{\link[grid]{grid-package}}.
#'
#' \bold{ggplot2-based plot:} A ggplot2 object is returned. Colors, axis limits
#' and other visual aspects can be added using standard ggplot2 functions
#' (\code{labs}, \code{scale_fill_continuous}, \code{labs}, etc.).
#'
#' @family plots
#'
#' @param x a matrix or an object of class \code{dist}.
#' @param order an object of class \code{ser_permutation} or the name of a seriation method.
#'   If \code{NULL} the order in \code{x} is plotted.
#' @param col a list of colors used. If \code{NULL}, a gray scale is used (for
#' matrix larger values are displayed darker and for \code{dist} smaller
#' distances are darker). For matrices containing logical data, black and white
#' is used. For matrices containing negative values a symmetric diverging color
#' palette is used.
#' @param main plot title.
#' @param xlab,ylab labels for the x and y axes.
#' @param zlim vector with two elements giving the range (min, max) for
#' representing the values in the matrix.
#' @param key logical; add a color key? No key is available for logical
#' matrices.
#' @param keylab string plotted next to the color key.
#' @param symkey logical; if \code{x} contains negative values, should the
#' color palate be symmetric (zero is in the middle)>
#' @param upper_tri,lower_tri,diag a logical indicating whether to show the
#' upper triangle, the lower triangle or the diagonal of the (distance) matrix.
#' @param row_labels,col_labels a logical indicating if row and column labels
#' in \code{x} should be displayed.  If \code{NULL} then labels are displayed
#' if the \code{x} contains the appropriate dimname and the number of labels is
#' 25 or less. A character vector of the appropriate length with labels can
#' also be supplied.
#' @param prop logical; change the aspect ratio so cells in the image have a
#' equal width and height.
#' @param flip_axes logical; exchange rows and columns for plotting.
#' @param reverse_columns logical; revers the order of how the columns are
#' displayed.
#' @param \dots if `order` is the name of a seriation method then further arguments are  passed
#'   on to the seriatation method, otherwise they are ignored.
#' @param newpage,pop,gp Start plot on a new page, pop the viewports after
#' plotting, and use the supplied \code{gpar} object (see \pkg{grid}).
#' @returns Nothing.
#'
#' @author Christian Buchta and Michael Hahsler
#' @keywords hplot
#' @examples
#' set.seed(1234)
#'
#' ## Example: Logical Matrix
#' x <- matrix(sample(c(FALSE, TRUE), 300, rep = TRUE), ncol = 10,
#'   dimnames = list(1:30, LETTERS[1:10]))
#'
#' # Matrix for logical values. TRUE values are dark and no color key is shown. There are too many
#' # Row labels (>25) so they are suppressed.
#' pimage(x)
#'
#' # Show all labels and flip axes or reverse columns
#' pimage(x, row_labels = TRUE, col_labels = TRUE, flip_axes = TRUE)
#' pimage(x, row_labels = TRUE, col_labels = TRUE, reverse_columns = TRUE)
#'
#' # Reorder matrix, use custom colors, and add a title.
#' pimage(x, order = seriate(x), row_labels = TRUE, col_labels = TRUE,
#'   col = c("white", "red"), main = "Random Data (Reordered)")
#'
#' ## Example: Positive Matrix
#' x <- matrix(runif(100), ncol = 10,
#'   dimnames = list(LETTERS[1:10], paste0("X", 1:10)))
#'
#' pimage(x, prop = TRUE)
#'
#' ## Example: Pos/Neg. Matrix
#' x <- matrix(rnorm(100), ncol = 10,
#'   dimnames = list(LETTERS[1:10], paste0("X", 1:10)))
#'
#' pimage(x, prop = TRUE)
#'
#' ## Example: Distance Matrix
#' # Show a reordered distance matrix (distances between rows).
#' # Dark means low distance. The aspect ratio is automatically fixed to 1:1.
#' d <- dist(x)
#' pimage(d,  order = seriate(d),
#'   main = "Random Data (Distances)")
#'
#' # Supress the upper triangle and diagonal
#' pimage(d,  order = seriate(d), upper = FALSE, diag = FALSE,
#'   main = "Random Data (Distances)")
#'
#' # Show only distances that are smaller than 4 using limits on z.
#' pimage(d, order = seriate(d),
#'   main = "Random Data (Distances + Theshold)", zlim = c(0, 4))
#'
#' ## Example: Correlation Matrix
#' # we calculate correlation between rows and seriate the matrix
#' r <- cor(t(x))
#' r <- permute(r, seriate(r))
#' pimage(r, upper = FALSE, diag = FALSE, zlim = c(-1, 1), reverse_columns = TRUE,
#'   main = "Random Data (Correlation)")
#'
#' # Add to the plot using functions in package grid
#' # Note: pop = FALSE allows us to manipulate viewports
#' library("grid")
#' pimage(x, pop = FALSE)
#'
#' # available viewports are: "main", "colorkey", "plot", "image"
#' current.vpTree()
#'
#' # Highlight cell column 7 (G) / row 5 (from top)/col with a red arrow starting at 5/2
#' # Note: columns are x and rows are y.
#' downViewport(name = "image")
#' grid.lines(x = c(5, 7), y = c(2, 5), arrow = arrow(),
#'   default.units = "native", gp = gpar(col = "red", lwd = 3))
#'
#' # add a red box around rows 15 and 16
#' grid.rect(x = 0.5, y = 5.5, width = ncol(x), height = 2,
#'   just = "left",
#'   default.units = "native", gp = gpar(col = "red", lwd = 3, fill = NA))
#'
#' ## remove the viewports
#' popViewport(0)
#'
#' ## put several pimages on a page (use grid viewports and newpage = FALSE)
#' # set up grid layout
#' library(grid)
#' grid.newpage()
#' top_vp <- viewport(layout = grid.layout(nrow = 1, ncol = 2,
#'   widths = unit(c(.4, .6), unit = "npc")))
#' col1_vp <- viewport(layout.pos.row = 1, layout.pos.col = 1, name = "col1_vp")
#' col2_vp <- viewport(layout.pos.row = 1, layout.pos.col = 2, name = "col2_vp")
#' splot <- vpTree(top_vp, vpList(col1_vp, col2_vp))
#' pushViewport(splot)
#'
#' seekViewport("col1_vp")
#' o <- seriate(x)
#' pimage(x, o, col_labels = FALSE, main = "Random Data",
#'   newpage = FALSE)
#'
#' seekViewport("col2_vp")
#' ## add the reordered dissimilarity matrix for rows
#' d <- dist(x)
#' pimage(d, o[[1]], main = "Random Data",
#'   newpage = FALSE)
#'
#' popViewport(0)
#'
#' ##-------------------------------------------------------------
#' ## ggplot2 Examples
#' if (require("ggplot2")) {
#'
#' library("ggplot2")
#'
#' ## Example: Logical Matrix
#' x <- matrix(sample(c(FALSE, TRUE), 300, rep = TRUE), ncol = 10,
#'   dimnames = list(1:30, LETTERS[1:10]))
#'
#' # Matrix for logical values. TRUE values are dark. There are too many
#' # Row labels (>25) so they are suppressed.
#' ggpimage(x)
#'
#' # Show all labels and flip axes or reverse columns
#' ggpimage(x, flip_axes = TRUE, row_labels = TRUE, col_labels = TRUE)
#' ggpimage(x, reverse_columns = TRUE, row_labels = TRUE, col_labels = TRUE)
#'
#' # Add lines
#' ggpimage(x) +
#'   geom_hline(yintercept = seq(0, nrow(x)) + .5) +
#'   geom_vline(xintercept = seq(0, ncol(x)) + .5)
#'
#' # Reorder matrix, use custom colors, add a title,
#' # and hide colorkey.
#' ggpimage(x, order = seriate(x), row_labels = TRUE, col_labels = TRUE) +
#'   scale_fill_manual(values = c("grey90", "red")) +
#'   theme(legend.position = "none") +
#'   labs(title = "Random Data")
#'
#' ## Example: Positive Matrix
#' x <- matrix(runif(100), ncol = 10,
#'   dimnames = list(LETTERS[1:10], paste0("X", 1:10)))
#'
#' ggpimage(x, order = seriate(x)) +
#'   labs(title = "Random Data")
#'
#' #' ## Example: Pos/Neg. Matrix
#' x <- matrix(rnorm(100), ncol = 10,
#'   dimnames = list(LETTERS[1:10], paste0("X", 1:10)))
#'
#' ggpimage(x, order = seriate(x)) +
#'   labs(title = "Random Data")
#'
#' ## Example: Distance Matrix
#' # Show a reordered distance matrix (distances between rows).
#' # Dark means low distance. The aspect ratio is automatically fixed to 1:1.
#' # The upper triangle is suppressed triangle
#' d <- dist(x)
#' ggpimage(d, order = seriate(d)) +
#'   labs(title = "Random Data", subtitle = "Distances")
#'
#' # Hide upper triangle and diagonal
#' ggpimage(d, order = seriate(d), upper_tri = FALSE, diag = FALSE) +
#'   labs(title = "Random Data", subtitle = "Distances")
#'
#' # Show only distances that are smaller than 4 using limits on fill.
#' ggpimage(d,  order = seriate(d), zlim = c(0, 4)) +
#'   labs(title = "Random Data (Distances + Theshold)")
#'
#' ## Example: Correlation Matrix
#' # we calculate correlation between rows and seriate the matrix
#' r <- cor(t(x))
#' r <- permute(r, seriate(r))
#' ggpimage(r,  zlim = c(-1, 1), upper = FALSE, diag = FALSE, reverse_columns = TRUE) +
#'   geom_text(aes(x = col, y = row, label = round(x, 2)), color = "black", size = 4) +
#'   labs(title = "Random Data", subtitle = "Correlation")
#'
#' ## Example: Custom themes and colors
#' # Use ggplot2 themes with theme_set
#' old_theme <- theme_set(theme_linedraw())
#' ggpimage(d, order = seriate(d)) +
#'   labs(title = "Random Data (Distances)")
#' theme_set(old_theme)
#'
#' # Use custom color palettes: Gray scale, Colorbrewer (provided in ggplot2) and colorspace
#' ggpimage(d, order = seriate(d), upper_tri = FALSE) +
#'   scale_fill_gradient(low = "black", high = "white", na.value = "white")
#'
#' ggpimage(d, order = seriate(d), upper_tri = FALSE) +
#'   scale_fill_distiller(palette = "Spectral", direction = +1, na.value = "white")
#'
#' ggpimage(d, order = seriate(d), upper_tri = FALSE) +
#'   colorspace::scale_fill_continuous_sequential("Reds", rev = FALSE, na.value = "white")
#' }
#' @export
pimage <-
  function(x,
    order = NULL,
    ...)
UseMethod("pimage")

### Note for matrix large values are dark, for dist large values are light!
#' @rdname pimage
#' @export
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
    diag = TRUE,
    row_labels = NULL,
    col_labels = NULL,
    prop = FALSE,
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
      x <- permute(x, order, ...)

    # mask triangles
    if (any(!upper_tri ||
        !lower_tri ||
        !diag) &&
        nrow(x) != ncol(x))
      stop("Upper triange, lower triangle or diagonal can only be suppressed for square matrices!")
    if (!upper_tri)
      x[upper.tri(x)] <- NA
    if (!lower_tri)
      x[lower.tri(x)] <- NA
    if (!diag)
      diag(x) <- NA

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
       down <- downViewport("colorkey")
      .grid_colorkey(zlim,
        col = col,
        horizontal = FALSE,
        lab = keylab)
      upViewport(down)

    } else
      .grid_basic_layout(
        main = main,
        left = left_mar,
        right = unit(0, "lines"),
        bottom = bottom_mar,
        gp = gp
      )

    down <- downViewport("plot")
    .grid_image(
      x,
      col = col,
      zlim = zlim,
      prop = prop
    ) #, gp=gp)
    upViewport(down)

    ## axes and labs
    down <- downViewport("image")
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

    # it is always 2 up from main
    seekViewport("main")
    down <- 2

    if (pop)
      popViewport(down)
    else
      upViewport(down)
  }

#' @export
pimage.default <- pimage.matrix

# as.matrix does not work for table!
table2matrix <- function(M)
  matrix(M, ncol = ncol(M), dimnames = dimnames(M))

#' @rdname pimage
#' @export
pimage.table <- function(x, order = NULL, ...)
  pimage.matrix(table2matrix(x), order = order, ...)

#' @rdname pimage
#' @export
pimage.data.frame <- function(x, order = NULL, ...)
  pimage.matrix(as.matrix(x), order = order, ...)


## small values are dark
#' @rdname pimage
#' @export
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
    upper_tri = TRUE,
    lower_tri = TRUE,
    diag = TRUE,
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
      x <- permute(x, order, ...)

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
      diag = diag,
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
