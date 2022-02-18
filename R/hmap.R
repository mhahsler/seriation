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

#' Plot Heat Map Reordered Using Seriation
#'
#' Provides heatmaps reordered using several different seriation methods. This
#' includes dendrogram based reordering with optimal leaf order and matrix
#' seriation-based heat maps.
#'
#' For dendrogram based heat maps, the arguments are passed on to
#' [stats::heatmap()] in \pkg{stats}. The following arguments for `heatmap()`
#' cannot be used: `margins`, `Rowv`, `Colv`, `hclustfun`, `reorderfun`.
#'
#' For seriation-based heat maps further arguments include:
#' - `gp` an object of class `gpar` containing graphical
#'   parameters (see [gpar()] in package \pkg{grid}).
#' - `newpage` a logical indicating whether to start plot on a new
#'   page
#' - `prop` a logical indicating whether the height and width of `x` should
#'   be plotted proportional to its dimensions.
#' - `showdist` Display seriated dissimilarity matrices? Values are
#'   `"none"`, `"both"`, `"rows"` or `"columns"`.
#' - `key` logical; show a colorkey?
#' - `key.lab` Label plotted next to the color key.
#' - `margins` bottom and right-hand-side margins are calculated
#'   automatically or can be specifies as a vector of two numbers (in lines).
#' - `zlim` range of values displayed.
#' - `col`, `col_dist` color palettes used.
#'
#' @family plots
#'
#' @param x a matrix or a dissimilarity matrix of class dist. If a
#' dissimilarity matrix is used, then the `distfun` is ignored.
#' @param distfun function used to compute the distance (dissimilarity) between
#' both rows and columns. For `gghmap()`, this
#' parameter is passed on in `control`.
#' @param method a character strings indicating the used seriation algorithm
#' (see [seriate.dist()]). If the method results in a dendrogram then
#' [stats::heatmap()] is used to show the dendrograms, otherwise
#' reordered distance matrices are shown instead.
#' @param control a list of control options passed on to the seriation
#' algorithm specified in `method`.
#' @param scale character indicating if the values should be centered and
#' scaled in either the row direction or the column direction, or none. Default
#' is none.
#' @param showDend Show dendrograms in the margin?
#' @param col a list of colors used.
#' @param row_labels,col_labels a logical indicating if row and column labels
#' in `x` should be displayed.  If `NULL` then labels are displayed
#' if the `x` contains the appropriate dimname and the number of labels is
#' 25 or less. A character vector of the appropriate length with labels can
#' also be supplied.
#' @param prop logical; change the aspect ratio so cells in the image have a
#' equal width and height.
#' @param \dots further arguments passed on to [stats::heatmap()].
#' @return An invisible list with elements:
#' \item{rowInd, colInd}{index permutation vectors.}
#' \item{reorder_method}{name of the method used to reorder the matrix.}
#'
#' The list may contain additional elements (dendrograms, colors, etc).
#'
#' @author Michael Hahsler
#' @keywords hplot
#' @examples
#' data("Wood")
#'
#' # regular heatmap from package stats
#' heatmap(Wood, main = "Wood (standard heatmap)")
#'
#' # Default heatmap does Euclidean distance, hierarchical clustering with
#' # average-link and optimal leaf ordering
#' hmap(Wood, main = "Wood (opt. leaf ordering)")
#'
#' # Heatmap shown without dendrograms (used pimage)
#' hmap(Wood, main = "Wood (opt. leaf ordering)",  showDend = FALSE)
#'
#' # Heatmap with correlation-based distance, green-red color (greenred is
#' # predefined) and optimal leaf ordering and no row label
#' dist_cor <- function(x) as.dist(1 - cor(t(x)))
#' hmap(Wood, distfun = dist_cor, col = greenred(100), row_labels = FALSE,
#'   main = "Wood (reorded by corr. between obs.)")
#'
#' # Heatmap with order based on the angle in two-dimensional MDS space.
#' hmap(Wood, method = "MDS_angle", col = greenred(100), row_labels = FALSE,
#'   main = "Wood (reorderd using ange in MDS space)")
#'
#' # Heatmap for distances
#' d <- dist(Wood)
#' hmap(d, method = "OLO", main = "Wood (Euclidean distances)")
#'
#' # order-based with dissimilarity matrices
#' hmap(Wood, method = "MDS_angle", showdist = "both",
#'   col = greenred(100), col_dist = greens(100),
#'   keylab = "norm. Expression", main = "Wood (reporderd with distances)")
#'
#' # Manually seriate and plot as pimage.
#' o <- seriate(Wood, method = "heatmap", control = list(dist_fun = dist, seriation_method = "OLO"))
#' o
#'
#' pimage(Wood, o, prop = FALSE)
#'
#' # Note: method heatmap calculates reorderd hclust objects which can be used for many heatmap
#' # implementations.
#' heatmap(Wood, Rowv = as.dendrogram(o[[1]]), Colv = as.dendrogram(o[[2]]))
#'
#' # ggplot 2 version does not support dendrograms in the margin (for now)
#' if (require("ggplot2")) {
#'   library("ggplot2")
#'
#'   gghmap(Wood) + labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   # More parameters (see ? ggpimage): reverse column order and flip axes, make a proportional plot
#'   gghmap(Wood, reverse_columns = TRUE) +
#'     labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   gghmap(Wood, flip_axes = TRUE) +
#'     labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   gghmap(Wood, flip_axes = TRUE, prop = TRUE) +
#'     labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   gghmap(Wood, distfun = dist_cor) +
#'     labs(title = "Wood", subtitle = "Reorded by correlation between observations") +
#'     scale_fill_gradient2(low = "darkgreen", high = "red")
#'
#'   gghmap(d, prop = TRUE) +
#'     labs(title = "Wood", subtitle = "Euclidean distances, reordered)")
#'
#'   # Note: the ggplot2-based version cannot show distance matrices in the same plot.
#'
#'   # Manually seriate and plot as pimage.
#'   o <- seriate(Wood, method = "heatmap", control = list(dist_fun = dist,
#'     seriation_method = "OLO"))
#'   o
#'
#'   ggpimage(Wood, o, prop = FALSE)
#' }
#' @export
hmap <- function(x,
  distfun = stats::dist,
  method = "OLO",
  control = NULL,
  scale = c("none", "row", "column"),
  showDend = TRUE,
  col = NULL,
  row_labels = NULL,
  col_labels = NULL,
  ...) {
  scale <- match.arg(scale)

  if (is.null(col)) {
    if (any(x < 0, na.rm = TRUE))
      col <- .diverge_pal()
    else
      col <- .sequential_pal()
  }

  # dist or matrix?
  if (inherits(x, "dist")) {
    dist_row <- dist_col <- x
    o_col <- o_row <- seriate(x,
      method = method, control = control)[[1]]
    x <- as.matrix(x)

    # dist uses reversed colors!
    col <- rev(col)
  } else {
    if (!is.matrix(x))
      x <- as.matrix(x)

    if (!is.null(scale)) {
      if (scale == "row")
        x <- t(scale(t(x)))
      if (scale == "col")
        x <- scale(x)
    }

    dist_row <- distfun(x)
    o_row <- seriate(dist_row,
      method = method, control = control)[[1]]

    #o_row <- ser_align(list(ser_permutation_vector(order(rowMeans(x, na.rm = TRUE), decreasing = TRUE)), o_row))[[2]]

    dist_col <- distfun(t(x))
    o_col <- seriate(dist_col,
      method = method, control = control)[[1]]

    #o_col <- ser_align(list(ser_permutation_vector(order(colMeans(x, na.rm = TRUE), decreasing = FALSE)), o_col))[[2]]
  }


  # is hierarchical? then let's do a heatmap from stats
  if (inherits(o_col, "hclust") && showDend) {
    # heatmap by default scales rows: we don't want that!
    # options are ignored for now: we use ...

    stats::heatmap(
      x,
      Rowv = stats::as.dendrogram(o_row),
      Colv = stats::as.dendrogram(o_col),
      scale = "none",
      col = col,
      labRow = row_labels,
      labCol = col_labels,
      ...
    )

  } else {
    ### we plot seriated distance matrices
    .hmap_dist(x, method, dist_row, dist_col, o_row, o_col, col = col,
      row_labels = row_labels, col_labels = col_labels, ...)
  }

  ## return permutation indices
  return(invisible(list(
    rowInd = o_row,
    colInd = o_col,
    seriation_method = method
  )))

}

## grid-based dissimilarity plot with seriation
.hmap_dist <-
  function(x,
    method,
    dist_row,
    dist_col,
    o_row,
    o_col,
    ...) {

    ## options
    options <- list(...)
    options <- .get_parameters(
      options,
      list(
        col       = if (any(x < 0))
          .diverge_pal()
        else
          .sequential_pal(),
        col_dist  = .sequential_pal(),
        prop      = FALSE,
        main      = NULL,
        key       = TRUE,
        keylab   = "",
        row_labels    = NULL,
        col_labels    = NULL,
        showdist  = "none",
        symm      = FALSE,
        margins   = NULL,
        zlim      = if (any(x < 0, na.rm = TRUE))
          max(abs(range(x, na.rm = TRUE))) * c(-1, 1)
        else
          range(x, na.rm = TRUE),
        newpage   = TRUE,
        gp        = gpar()
      )
    )

    options$col_dist <- rev(options$col_dist)
    .showdist_options <- c("none", "row", "column", "both")
    options$showdist <- .showdist_options[pmatch(options$showdist,
      .showdist_options)]
    if (is.na(options$showdist))
      stop("Unknown value for showdist. Use one of: ",
        paste(dQuote(.showdist_options), collapse = ", "))

    ## if symmetric then we only use o_row and dist_row
    if (length(o_row) == length(o_col) && options$symm == TRUE) {
      o_col <- o_row
      dist_col <- dist_row
    }

    x <- permute(x, ser_permutation(o_row, o_col))

    if (options$showdist == "none") {
      pimage(
        x,
        col = options$col,
        main = options$main,
        zlim = options$zlim,
        row_labels = options$row_labels,
        col_labels = options$col_labels,
        prop = options$prop,
        key = options$key,
        newpage = options$newpage,
        gp = options$gp
      )
      return()
    }

    dist_row <- permute(dist_row, o_row)
    dist_col <- permute(dist_col, o_col)

    # deal with row/col labels
    row_labels <- options$row_labels
    col_labels <- options$col_labels
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

    ## Note: we need a list to store units!
    if (is.null(options$margins)) {
      options$margins <- list(unit(1, "lines"), unit(1, "lines"))
      if (col_labels)
        options$margins[[1]] <-
          max(stringWidth(colnames(x))) + unit(2, "lines")
      if (row_labels)
        options$margins[[2]] <-
          max(stringWidth(rownames(x))) + unit(2, "lines")
      all_names <-
        c("", if (col_labels)
          colnames(x), if (row_labels)
            rownames(x))
      options$margins[[3]] <-
        max(stringWidth(all_names)) + unit(2, "lines")
    } else
      options$margins <-
      list(
        unit(options$margins[1], "lines"),
        unit(options$margins[2], "lines"),
        unit(max(options$margins), "lines")
      )


    ## plot
    if (options$newpage)
      grid.newpage()

    ## surrounding viewport
    pushViewport(viewport(
      layout = grid.layout(
        nrow = 3 ,
        ncol = 3,
        widths = unit.c(
          unit(1, "lines"),
          unit(1, "snpc") - options$margins[[3]] - unit(3, "lines"),
          options$margins[[2]]
        ),
        heights = unit.c(
          unit(3, "lines"),
          # main
          unit(1, "snpc") - options$margins[[3]] - unit(3, "lines"),
          options$margins[[1]]
        )
      ),
      width = unit(1, "snpc"),
      height = unit(1, "snpc"),
      gp = options$gp
    ))


    ## main title
    if (!is.null(options$main)) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
      grid.text(options$main, gp = gpar(cex = 1.3))
      upViewport(1)
    }


    ## plots
    if (options$prop) {
      widths <- unit.c(
        unit(1 - ncol(x) / sum(ncol(x), nrow(x)), "npc") - unit(.25, "lines"),
        unit(.5, "lines"),
        unit(ncol(x) / sum(ncol(x), nrow(x)), "npc") - unit(.25, "lines")
      )

      heights <- unit.c(
        unit(1 - nrow(x) / sum(ncol(x), nrow(x)), "npc") - unit(.25, "lines"),
        unit(.5, "lines"),
        #space
        unit(nrow(x) / sum(ncol(x), nrow(x)), "npc") - unit(.25, "lines")
      )
    } else{
      heights <- widths <- unit.c(unit(1, "null"),
        unit(.5, "lines"),   # space
        unit(1, "null"))
    }

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 3,
          ncol = 3,
          widths = widths,
          heights = heights
        ),
        width = unit(1, "snpc"),
        height = unit(1, "snpc"),
        layout.pos.row = 2,
        layout.pos.col = 2
      )
    )

    # data
    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3))

    .grid_image(x,
      col = options$col,
      gp = options$gp,
      zlim = options$zlim)

    downViewport("image")
    if (col_labels)
      grid.text(
        colnames(x),
        y = unit(-1, "lines"),
        x = unit(1:ncol(x), "native"),
        rot = 90,
        just = "right"
      ) # , gp=options$gp)
    if (row_labels)
      grid.text(
        rownames(x),
        x = unit(1, "npc") + unit(1, "lines"),
        y = unit(1:nrow(x), "native"),
        just = "left"
      ) #, gp=options$gp)
    popViewport(1)

    popViewport(1)

    # rows
    if (options$showdist %in% c("row", "both")) {
      pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
      .grid_image(as.matrix(dist_row),
        col = options$col_dist,
        gp = options$gp)
      popViewport(1)
    }

    # cols
    if (options$showdist %in% c("column", "both")) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
      .grid_image(as.matrix(dist_col),
        col = options$col_dist,
        gp = options$gp)
      popViewport(1)
    }

    # colorkey
    if (options$key) {
      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))

      pushViewport(viewport(
        width = unit(0.5, "npc"),
        height = unit(1, "lines")
      ))
      .grid_colorkey(
        options$zlim,
        col = options$col,
        lab = options$keylab,
        gp = options$gp
      )
      popViewport(2)
    }

    popViewport(2)

  }
