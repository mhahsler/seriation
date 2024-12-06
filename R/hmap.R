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
#' (see [seriate.dist()]).
#' If the method results in a dendrogram then
#' [stats::heatmap()] is used to show the dendrograms, otherwise
#' reordered distance matrices are shown instead.
#' @param control a list of control options passed on to the seriation
#' algorithm specified in `method`.
#' @param scale character indicating if the values should be centered and
#' scaled in either the row direction or the column direction, or none. Default
#' is none.
#' @param plot_margins character indicating what to show in the margins. Options are:
#'   `"auto"`, `"dendrogram"`, `"distances"`, or `"none"`.
#' @param col a list of colors used.
#' @param col_dist colors used for displaying distances.
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
#' # Default heatmap does Euclidean distance, hierarchical clustering with
#' # complete-link and optimal leaf ordering. Note that the rows are
#' # ordered top-down in the seriation order (stats::heatmap orders in reverse)
#' hmap(Wood, main = "Wood (opt. leaf ordering)")
#' hmap(Wood, plot_margins = "distances", main = "Wood (opt. leaf ordering)")
#' hmap(Wood, plot_margins = "none", main = "Wood (opt. leaf ordering)")
#'
#' # Heatmap with correlation-based distance, green-red color (greenred is
#' # predefined) and optimal leaf ordering and no row label
#' dist_cor <- function(x) as.dist(sqrt(1 - cor(t(x))))
#' hmap(Wood, distfun = dist_cor, col = greenred(100),
#'   main = "Wood (reorded by corr. between obs.)")
#'
#' # Heatmap for distances
#' d <- dist(Wood)
#' hmap(d, main = "Wood (Euclidean distances)")
#'
#' # order-based with dissimilarity matrices
#' hmap(Wood, method = "MDS_angle",
#'   col = greenred(100), col_dist = greens(100, power = 2),
#'   keylab = "norm. Expression", main = "Wood (reorderd with distances)")
#'
#' # Manually create a simple heatmap with pimage.
#' o <- seriate(Wood, method = "heatmap",
#'    control = list(dist_fun = dist, seriation_method = "OLO_ward"))
#' o
#'
#' pimage(Wood, o)
#'
#' # Note: method heatmap calculates reorderd hclust objects which can be used
#' #       for many heatmap implementations like the standard implementation in
#' #       package stats.
#' heatmap(Wood, Rowv = as.dendrogram(o[[1]]), Colv = as.dendrogram(o[[2]]))
#'
#' # ggplot 2 version does not support dendrograms in the margin (for now)
#' if (require("ggplot2")) {
#'   library("ggplot2")
#'
#'   gghmap(Wood) + labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   gghmap(Wood, flip_axes = TRUE, prop = TRUE) +
#'     labs(title = "Wood", subtitle = "Optimal leaf ordering")
#'
#'   dist_cor <- function(x) as.dist(sqrt(1 - cor(t(x))))
#'   gghmap(Wood, distfun = dist_cor) +
#'     labs(title = "Wood", subtitle = "Reorded by correlation between observations") +
#'     scale_fill_gradient2(low = "darkgreen", high = "red")
#'
#'   gghmap(d, prop = TRUE) +
#'     labs(title = "Wood", subtitle = "Euclidean distances, reordered")
#'
#'   # Note: the ggplot2-based version currently cannot show distance matrices
#'   #      in the same plot.
#'
#'   # Manually seriate and plot as pimage.
#'   o <- seriate(Wood, method = "heatmap", control = list(dist_fun = dist,
#'     seriation_method = "OLO_ward"))
#'   o
#'
#'   ggpimage(Wood, o)
#' }
#' @export
hmap <- function(x,
                 distfun = stats::dist,
                 method = "OLO_complete",
                 control = NULL,
                 scale = c("none", "row", "column"),
                 plot_margins = "auto",
                 col = NULL,
                 col_dist = grays(power = 2),
                 row_labels = NULL,
                 col_labels = NULL,
                 ...) {
  scale <- match.arg(scale)
  plot_margins <-
    match.arg(plot_margins, c("auto", "dendrogram", "distances", "none"))

  if (is.null(col)) {
    if (any(x < 0, na.rm = TRUE))
      col <- .diverge_pal()
    else
      col <- .sequential_pal()
  }

  # dist or matrix?
  if (inherits(x, "dist")) {
    dist_row <- dist_col <- x
    o <- seriate(x,
                 method = method, control = control)[[1]]
    o <- ser_permutation(o, o)
    x <- as.matrix(x)

    # dist uses reversed colors!
    col <- rev(col)
  } else {
    if (!is.matrix(x))
      x <- as.matrix(x)

    o <-
      seriate(
        x,
        "Heatmap",
        seriation_method = method,
        dist_fun = distfun,
        seriation_control = control,
        scale = scale
      )
  }

  if (plot_margins == "auto") {
    if (all(sapply(o, inherits, "hclust")))
      plot_margins <- "dendrogram"
    else
      plot_margins <- "distances"
  }

  if (plot_margins == "dendrogram" &&
      !all(sapply(o, inherits, "hclust"))) {
    warning(
      "Dendrogramms not available for all dimensions! Plotting distance matrices instead."
    )
    plot_margins <- "distances"
  }

  if (plot_margins == "dendrogram") {
    # heatmap by default scales rows: we don't want that!
    # options are ignored for now: we use ...

    stats::heatmap(
      x,
      Rowv = stats::as.dendrogram(rev(o[[1]])),
      Colv = stats::as.dendrogram(o[[2]]),
      scale = scale,
      col = col,
      labRow = row_labels,
      labCol = col_labels,
      ...
    )

  } else if (plot_margins == "distances") {
    ### we plot seriated distance matrices
    #pimage(x, o, col = col, row_labels = row_labels, col_labels = col_labels, ...)
    .hmap_dist(
      x,
      method,
      dist_row = distfun(x),
      dist_col = distfun(t(x)),
      o,
      col = col,
      col_dist = col_dist,
      row_labels = row_labels,
      col_labels = col_labels,
      ...
    )
  } else
    pimage(x,
           o,
           col = col,
           row_labels = row_labels,
           col_labels = col_labels,
           ...)

  ## return permutation indices
  return(invisible(list(
    o = o,
    seriation_method = method
  )))

}

## grid-based dissimilarity plot with seriation
.hmap_dist <-
  function(x,
           method,
           dist_row,
           dist_col,
           o,
           ...) {
    o_row <- o[[1]]
    o_col <- o[[2]]

    ## options
    options <- list(...)
    options <- .get_parameters(
      options,
      list(
        col       = if (any(x < 0))
          .diverge_pal()
        else
          .sequential_pal(),
        col_dist  = grays,
        prop      = FALSE,
        main      = NULL,
        key       = TRUE,
        keylab   = "",
        row_labels    = NULL,
        col_labels    = NULL,
        showdist  = "both",
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
