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

#' Dissimilarity Plot
#'
#' Visualizes a dissimilarity matrix using seriation and matrix shading using
#' the method developed by Hahsler and Hornik (2011). Entries with lower
#' dissimilarities (higher similarity) are plotted darker. Dissimilarity plots
#' can be used to uncover hidden structure in the data and judge cluster
#' quality.
#'
#' The plot can also be used to visualize cluster quality (see Ling 1973).
#' Objects belonging to the same cluster are displayed in consecutive order.
#' The placement of clusters and the within cluster order is obtained by a
#' seriation algorithm which tries to place large similarities/small
#' dissimilarities close to the diagonal. Compact clusters are visible as dark
#' squares (low dissimilarity) on the diagonal of the plot. Additionally, a
#' Silhouette plot (Rousseeuw 1987) is added. This visualization is similar to
#' CLUSION (see Strehl and Ghosh 2002), however, allows for using arbitrary
#' seriating algorithms.
#'
#' **Note:** Since [pimage()] uses \pkg{grid}, it should not be mixed
#' with base R primitive plotting functions.
#'
#' @family plots
#'
#' @param x an object of class [dist].
#' @param labels `NULL` or an integer vector of the same length as
#' rows/columns in `x` indicating the cluster membership for each object
#' in `x` as consecutive integers starting with one. The labels are used
#' to reorder the matrix.
#' @param method A single character string indicating the seriation method used
#' to reorder the clusters (inter cluster seriation) as well as the objects
#' within each cluster (intra cluster seriation).  If different algorithms for
#' inter and intra cluster seriation are required, `method` can be a
#' `list` of two named elements (`inter_cluster` and
#' `intra_cluster` each containing the name of the respective seriation
#' method. Use [list_seriation_methods()] with `kind = "dist"` to find available algorithms.
#'
#' Set method to `NA` to plot the matrix as is (no or, if cluster labels
#' are supplied, only coarse seriation). For intra cluster reordering with the
#' special method `"silhouette width"` is available (for `dissplot()`
#' only). Objects in clusters are then ordered by silhouette width (from
#' silhouette plots).  If no `method` is given, the default method of
#' [seriate.dist()] is used.
#'
#' A third list element (named `aggregation`) can be added to control how
#' inter cluster dissimilarities are computed from from the given dissimilarity
#' matrix. The choices are `"avg"` (average pairwise dissimilarities;
#' average-link), `"min"` (minimal pairwise dissimilarities; single-link),
#' `"max"` (maximal pairwise dissimilarities; complete-link), and
#' `"Hausdorff"` (pairs up each point from one cluster with the most
#' similar point from the other cluster and then uses the largest dissimilarity
#' of paired up points).
#' @param control a list of control options passed on to the seriation
#' algorithm.  In case of two different seriation algorithms, `control`
#' can contain a list of two named elements (`inter_cluster` and
#' `intra_cluster`) containing each a list with the control options for
#' the respective algorithm.
#' @param upper_tri,lower_tri,diag a logical indicating whether to show the upper triangle, the
#' lower triangle or the diagonal of the distance matrix. The string "average" can also be used
#' to display within and between cluster averages in the two triangles.
#' @param cluster_labels a logical indicating whether to display cluster labels
#' in the plot.
#' @param cluster_lines a logical indicating whether to draw lines to separate
#' clusters.
#' @param reverse_columns a logical indicating if the clusters are displayed on
#' the diagonal from north-west to south-east (`FALSE`; default) or from
#' north-east to south-west (`TRUE`).
#' @param options a list with options for plotting the matrix (`dissplot`
#' only).
#' - `plot` a logical indicating if a plot should
#'   be produced.  if `FALSE`, the returned object can be plotted later
#'   using the function `plot` which takes as the second argument a list of
#'   plotting options (see `options` below).
#' - `silhouettes` a logical indicating whether to include a silhouette plot
#'   (see Rousseeuw, 1987).
#' - `threshold` a numeric. If used, only plot distances
#'   below the threshold are displayed. Consider also using `zlim` for this
#'   purpose.
#' - `col` colors used for the image plot.
#' - `key` a logical indicating whether to place a color key below the plot.
#' - `zlim` range of values to display (defaults to range `x`).
#' - `axes` `"auto"` (default; enabled for less than 25 objects), `"y"` or `"none"`.
#' - `main` title for the plot.
#' - `newpage` a logical indicating whether to start plot on a new page
#'   (see [grid.newpage()].
#' - `pop` a logical indicating whether to pop the created viewports?
#'   (see package \pkg{grid})
#' - `gp`, `gp_lines`, `gp_labels` objects of class `gpar` containing graphical parameters for the plot
#'   lines and labels (see [gpar()].
#' @param ... `dissplot()`: further arguments are added to `options`.
#'   `ggdissplot()` further arguments are passed on to [ggpimage()].
#' @return `dissplot()` returns an invisible object of class
#' `cluster_proximity_matrix` with the following elements:
#' \item{order}{`NULL` or integer vector giving the order used to plot `x`.}
#' \item{cluster_order}{ `NULL` or integer vector giving the order of the
#'   clusters as plotted.}
#' \item{method}{ vector of character strings indicating
#'   the seriation methods used for plotting `x`.}
#' \item{k}{ `NULL` or integer scalar giving the number of clusters generated.}
#' \item{description}{ a `data.frame` containing information (label, size, average
#'   intra-cluster dissimilarity and the average silhouette) for the clusters as
#'   displayed in the plot (from top/left to bottom/right).}
#'
#' This object can be used for plotting via `plot(x, options = NULL, ...)`,
#' where `x` is the object and `options` contains a list with
#' plotting options (see above).
#'
#' `ggdissplot()` returns a ggplot2 object representing the plot.
#'
#' @returns The plot description as an object of class `reordered_cluster_dissimilarity_matrix`.
#'
#' @author Michael Hahsler
#' @references
#' Hahsler, M. and Hornik, K. (2011): Dissimilarity plots: A visual
#' exploration tool for partitional clustering. \emph{Journal of Computational
#' and Graphical Statistics,} \bold{10}(2):335--354.
#' \doi{10.1198/jcgs.2010.09139}
#'
#' Ling, R.F. (1973): A computer generated aid for cluster analysis.
#' \emph{Communications of the ACM,} \bold{16}(6), 355--361.
#' \doi{10.1145/362248.362263}
#'
#' Rousseeuw, P.J. (1987): Silhouettes: A graphical aid to the interpretation
#' and validation of cluster analysis. \emph{Journal of Computational and
#' Applied Mathematics,} \bold{20}(1), 53--65.
#' \doi{10.1016/0377-0427(87)90125-7}
#'
#' Strehl, A. and Ghosh, J. (2003): Relationship-based clustering and
#' visualization for high-dimensional data mining. \emph{INFORMS Journal on
#' Computing,} \bold{15}(2), 208--230.
#' \doi{10.1287/ijoc.15.2.208.14448}
#' @keywords hplot cluster
#' @examples
#' data("iris")
#'
#' # shuffle rows
#' x_iris <- iris[sample(seq(nrow(iris))), -5]
#' d <- dist(x_iris)
#'
#' # Plot original matrix
#' dissplot(d, method = NA)
#'
#' # Plot reordered matrix using the nearest insertion algorithm (from tsp)
#' dissplot(d, method = "TSP", main = "Seriation (TSP)")
#'
#' # Cluster iris with k-means and 3 clusters and reorder the dissimality matrix
#' l <- kmeans(x_iris, centers = 3)$cluster
#' dissplot(d, labels = l, main = "k-means")
#'
#' # show only distances as lower triangle
#' dissplot(d, labels = l, main = "k-means", lower_tri = TRUE, upper_tri = FALSE)
#'
#' # Use a grid layout to place several plots on a page
#' library("grid")
#' grid.newpage()
#' pushViewport(viewport(layout=grid.layout(nrow = 2, ncol = 2),
#'     gp = gpar(fontsize = 8)))
#' pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
#'
#' # Visualize the clustering (using Spectral between clusters and MDS within)
#' res <- dissplot(d, l, method = list(inter = "Spectral", intra = "MDS"),
#'   main = "K-Means + Seriation", newpage = FALSE)
#'
#' popViewport()
#' pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
#'
#' # More visualization options. Note that we reuse the reordered object res!
#' # color: use 10 shades red-blue, biased towards small distances
#' plot(res, main = "K-Means + Seriation (red-blue + biased)",
#'     col= bluered(10, bias = .5), newpage = FALSE)
#'
#' popViewport()
#' pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
#'
#' # Threshold (using zlim) and cubic scale to highlight differences
#' plot(res, main = "K-Means + Seriation (cubic + threshold)",
#'     zlim = c(0, 2), col = grays(100, power = 3), newpage = FALSE)
#'
#' popViewport()
#' pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
#'
#' # Use gray scale with logistic transformation
#' plot(res, main = "K-Means + Seriation (logistic scale)",
#'   col = gray(
#'     plogis(seq(max(res$x_reordered), min(res$x_reordered), length.out = 100),
#'       location = 2, scale = 1/2, log = FALSE)
#'     ),
#'   newpage = FALSE)
#'
#' popViewport(2)
#'
#' # The reordered_cluster_dissimilarity_matrix object
#' res
#' names(res)
#'
#' # ggplot-based dissplot
#' if (require("ggplot2")) {
#'   library("ggplot2")
#'
#'   # Plot original matrix
#'   ggdissplot(d, method = NA)
#'
#'   # Plot seriated matrix
#'   ggdissplot(d, method = "TSP") +
#'     labs(title = "Seriation (TSP)")
#'
#'   # Cluster iris with k-means and 3 clusters
#'   l <- kmeans(x_iris, centers = 3)$cluster
#'
#'   ggdissplot(d, labels = l) +
#'     labs(title = "K-means + Seriation")
#'
#'   # show only lower triangle
#'   ggdissplot(d, labels = l, lower_tri = TRUE, upper_tri = FALSE) +
#'     labs(title = "K-means + Seriation")
#'
#'   # No lines or cluster labels and add a label for the color key (fill)
#'   ggdissplot(d, labels = l, cluster_lines = FALSE, cluster_labels = FALSE) +
#'     labs(title = "K-means + Seriation", fill = "Distances\n(Euclidean)")
#'
#'   # Diverging color palette with manual set midpoint and different seriation methods
#'   ggdissplot(d, l, method = list(inter = "Spectral", intra = "MDS")) +
#'     labs(title = "K-Means + Seriation", subtitle = "biased color scale") +
#'     scale_fill_gradient2(midpoint = median(d))
#'
#'   # Use manipulate scale using package scales
#'   library("scales")
#'
#'   # Threshold (using limit and na.value) and cubic scale to highlight differences
#'   cubic_dist_trans <- trans_new(
#'     name = "cubic",
#'     # note that we have to do the inverse transformation for distances
#'     trans = function(x) x^(1/3),
#'     inverse = function(x) x^3
#'   )
#'
#'   ggdissplot(d, l, method = list(inter = "Spectral", intra = "MDS")) +
#'     labs(title = "K-Means + Seriation", subtitle = "cubic + biased color scale") +
#'     scale_fill_gradient(low = "black", high = "white",
#'       limit = c(0,2), na.value = "white",
#'       trans = cubic_dist_trans)
#'
#'   # Use gray scale with logistic transformation
#'   logis_2_.5_dist_trans <- trans_new(
#'     name = "Logistic transform (location, scale)",
#'     # note that we have to do the inverse transformation for distances
#'     trans = function(x) plogis(x, location = 2, scale = .5, log = FALSE),
#'     inverse = function(x) qlogis(x, location = 2, scale = .5, log = FALSE),
#'   )
#'
#'   ggdissplot(d, l, method = list(inter = "Spectral", intra = "MDS")) +
#'     labs(title = "K-Means + Seriation", subtitle = "logistic color scale") +
#'     scale_fill_gradient(low = "black", high = "white",
#'       trans = logis_2_.5_dist_trans,
#'       breaks = c(0, 1, 2, 3, 4))
#' }
#' @export
dissplot <- function(x,
  labels = NULL,
  method = "Spectral",
  control = NULL,
  lower_tri = TRUE,
  upper_tri = "average",
  diag = TRUE,
  cluster_labels = TRUE,
  cluster_lines = TRUE,
  reverse_columns = FALSE,
  options = NULL,
  ...) {
  ## add ... to options
  options <- c(options, list(...))
  options$cluster_labels <- cluster_labels
  options$cluster_lines <- cluster_lines
  options$reverse_columns <- reverse_columns

  ## make x dist
  if (!inherits(x, "dist")) {
    if (is.matrix(x) && isSymmetric(x))
      x <- as.dist(x)
    else
      stop("Argument 'x' cannot safely be coerced to class 'dist'.")
  }

  a <- .arrange_dissimilarity_matrix(x,
    labels = labels,
    method = method,
    control = control)

  if (is.null(options$plot) || options$plot)
    plot(a, lower_tri, upper_tri, diag, options)

  invisible(a)
}


## work horse
.arrange_dissimilarity_matrix <-
  function(x,
    labels = NULL,
    method = NULL,
    control = NULL) {
    ## x is already of class dist
    dim <- attr(x, "Size")
    diss_measure <- attr(x, "method")

    ## check labels
    if (!is.null(labels) && length(labels) != dim)
      stop("Number of labels in 'labels' does not match dimensions of 'x'.")

    m <- method


    ## set everything to NULL first
    order               <- NULL
    k                   <- NULL             # number of clusters
    sil                 <- NULL
    avgSil              <- NULL
    labels_unique       <- NULL
    cluster_dissimilarities <- NULL
    ## method$a means method$ aggregation (default is avg)
    aggregation	<- "avg"
    if (class(method) == "list" &&
        !is.null(method$a))
      aggregation <- method$a

    if (class(method) != "list")
      method <-
      list(inter_cluster = m, intra_cluster = m)
    m <-
      pmatch(names(method),
        c("inter_cluster", "intra_cluster", "aggregation"))
    if (any(is.na(m)))
      stop("Unknown method component. Use 'inter_cluster', 'intra_cluster' and 'aggregation'.")
    names(method) <-
      c("inter_cluster", "intra_cluster", "aggregation")[m]

    if (class(control[[1]]) != "list") {
      control <- list(inter_cluster = control, intra_cluster = control)
    }

    if (!is.null(method$inter_cluster)
      && is.na(method$inter_cluster)) {
      ## no setiation
      if (!is.null(labels)) {
        ## do coarse seriation
        order <- order(labels)
        k <- length(unique(labels))
        ## calculate cluster_dissimilarities for later
        cluster_dissimilarities <- .cluster_dissimilarity(x, labels,
          aggregation)
        aggregation <- attr(cluster_dissimilarities, "method")
        ## calculate silhouette values for later use
        sil <- cluster::silhouette(labels, x)

      }
      ## else keep the matrix as is -- do not reorder

    } else if (is.null(labels)) {
      ## reorder whole matrix if no labels are given
      order <- seriate(x,
        method = method$inter_cluster,
        control = control$inter)[[1]]

      method$inter_cluster <- if (!is.null(attr(order, "method")))
        attr(order, "method")
      else
        method$inter_cluster

      order <- get_order(order)

    } else{
      ## reorder clusters for given labels
      ## get number of clusters k
      k <- length(unique(labels))

      ## reorder with average pairwise dissimilarites between clusters
      cluster_dissimilarities <- .cluster_dissimilarity(x, labels,
        aggregation)
      aggregation <- attr(cluster_dissimilarities, "method")

      if (k > 2) {
        cluster_order <- seriate(
          as.dist(cluster_dissimilarities),
          method = method$inter_cluster,
          control = control$inter
        )[[1]]

        method$inter_cluster <-
          if (!is.null(attr(cluster_order, "method")))
            attr(cluster_order, "method")
        else
          method$inter_cluster

        cluster_order <- get_order(cluster_order)
      } else{
        cluster_order <- 1:k
      }

      ## calculate silhouette values for later use
      sil <- cluster::silhouette(labels, x)

      ## determine order for matrix from cluster order
      order <- c()

      if (!is.null(method$intra_cluster) &&
          is.na(method$intra_cluster)) {
        ## no intra cluster ordering
        for (i in 1:k) {
          order <- c(order, which(labels == cluster_order[i]))
        }
        ##method$intra_cluster <- NA

      } else{
        ## intra cluster order

        for (i in 1:k) {
          take <- which(labels == cluster_order[i])

          ## only reorder for >1 elements
          if (length(take) > 1) {
            if (is.character(method$intra_cluster) &&
                match(
                  tolower(method$intra_cluster),
                  c("sil", "silhouette", "silhouette width"),
                  nomatch = 0
                ) > 0) {
              intra_order <-  order(sil[take, "sil_width"],
                decreasing = TRUE)

              method$intra_cluster <- "silhouette width"
            } else{
              ## we use .rearrange_dist instead of permute
              ## since we take only a subset!
              block <- .rearrange_dist(x, take)

              intra_order <-
                seriate(block,
                  method = method$intra_cluster,
                  control = control$intra)[[1]]

              method$intra_cluster <-
                if (!is.null(attr(intra_order, "method")))
                  attr(intra_order, "method")
              else
                method$intra_cluster

              intra_order <- get_order(intra_order)
            }

            order <- c(order, take[intra_order])

          } else{
            order <- c(order, take)
          }

        }
      }


      ## reorder cluster_dissimilarities for later
      cluster_dissimilarities  <-
        cluster_dissimilarities[cluster_order, cluster_order]

    }

    ## reorder matrix
    if (!is.null(order)) {
      x_reordered <- permute(x, order)
      labels <- labels[order]
    }
    else
      x_reordered <- x

    ## prepare for return value
    cluster_description <- NULL

    if (!is.null(labels)) {
      labels_unique   <-  unique(labels)

      ## reorder silhouettes
      sil <- sil[order,]

      ## calculate avg silhouettes
      avgSil <- sapply(labels_unique, function(x)
        mean(sil[sil[, "cluster"] == x, "sil_width"]))

      ## generate description
      cluster_description = data.frame(
        position        = c(1:k),
        label           = labels_unique,
        size            = tabulate(labels)[labels_unique],
        ## FIXME: this is not the average anymore!
        aggregated_dissimilarity =
          diag(cluster_dissimilarities)[labels_unique],
        avg_silhouette_width = avgSil
      )
    }

    ## clean order from names, etc.
    attributes(order) <- NULL

    structure(
      list(
        x_reordered        = x_reordered,
        labels             = labels,
        seriation_methods  = method,
        aggregation_method = aggregation,
        k                  = k,
        cluster_dissimilarities =  cluster_dissimilarities,
        sil                = sil,
        order              = order,
        cluster_order      = labels_unique,
        diss_measure       = diss_measure,
        description        =  cluster_description
      ),
      class = "reordered_cluster_dissimilarity_matrix"
    )
  }


## create panels with avg. dissimilarity
## a is an arrangement
.average_tri <- function(a,
  lower_tri = "average",
  upper_tri = TRUE,
  diag = TRUE) {
  if (!inherits(a, "reordered_cluster_dissimilarity_matrix"))
    stop("a needs to be a reordered_cluster_dissimilarity_matrix")

  upper_avg <- !is.na(pmatch(tolower(upper_tri), "average"))
  lower_avg <- !is.na(pmatch(tolower(lower_tri), "average"))

  k <- a$k
  labels <- a$labels
  labels_unique <- a$cluster_order
  cluster_dissimilarities <- a$cluster_dissimilarities
  m <- as.matrix(a$x_reordered)

  ## blank out if FALSE or NA
  if (is.na(upper_tri) || (is.logical(upper_tri) && !upper_tri)) {
    m[upper.tri(m)] <- NA
    upper_tri <- FALSE
  }
  if (is.na(lower_tri) || (is.logical(lower_tri) && !lower_tri)) {
    m[lower.tri(m)] <- NA
    lower_tri <- FALSE
  }

  ## do off-diagonal averages by cluster
  if (!is.null(cluster_dissimilarities) &&
      !is.null(labels) && (upper_avg || lower_avg)) {
    for (i in seq(2, k)) {
      for (j in seq(i - 1)) {
        ## check empty clusters
        if (is.na(labels_unique[i]))
          next
        if (is.na(labels_unique[j]))
          next

        ## lower panels
        if (lower_avg) {
          m[labels == labels_unique[i], labels == labels_unique[j]] <-
            cluster_dissimilarities[i, j]
        }

        ## upper panels
        if (upper_avg) {
          m[labels == labels_unique[j], labels == labels_unique[i]] <-
            cluster_dissimilarities[i, j]
        }
      }
    }


    ## do diagonal
    for (i in seq(1, k)) {
      block <- m[labels == labels_unique[i],
        labels == labels_unique[i]]

      if (upper_avg) {
        block[upper.tri(block, diag = TRUE)] <-
          cluster_dissimilarities[i, i]

        m[labels == labels_unique[i],
          labels == labels_unique[i]] <- block
      }

      if (lower_avg) {
        block[lower.tri(block, diag = TRUE)] <-
          cluster_dissimilarities[i, i]

        m[labels == labels_unique[i],
          labels == labels_unique[i]] <- block
      }
    }
  }

  if (!diag)
    diag(m) <- NA

  m
}


## plot for reordered_cluster_dissimilarity_matrix
#' @rdname dissplot
#' @export
plot.reordered_cluster_dissimilarity_matrix <-
  function(x,
    lower_tri = TRUE,
    upper_tri = "average",
    diag = TRUE,
    options = NULL,
    ...) {
    ## add ... to options
    options <- c(options, list(...))

    k       <- x$k
    dim     <- attr(x$x_reordered, "Size")
    labels  <- x$labels
    #labels_unique <- unique(labels)
    labels_unique <- x$cluster_order

    m  <- .average_tri(x,
      lower_tri = lower_tri,
      upper_tri = upper_tri,
      diag = diag)

    ## default plot options
    options <- .get_parameters(
      options,
      list(
        cluster_labels = TRUE,
        cluster_lines  = TRUE,
        reverse_columns	= FALSE,
        silhouettes    = FALSE,
        col            = NULL,
        threshold      = NULL,
        zlim           = NULL,
        key            = TRUE,
        main           = NULL,
        axes           = "auto",
        gp             = gpar(),
        gp_lines       = gpar(),
        gp_labels      = gpar(),
        newpage        = TRUE,
        pop            = TRUE
      )
    )

    if (is.null(options$col))
      options$col <- rev(.sequential_pal())
    else
      options$col <- rev(options$col)

    i <- pmatch(options$axes, c("auto", "x", "y", "both", "none"))
    if (is.na(i))
      stop("Illegal vaule for axes. Use: 'auto', 'x', 'y', 'both' or 'none'!")
    options$axes <- c("auto", "x", "y", "both", "none")[i]

    ## clear page
    if (options$newpage)
      grid.newpage()

    ## do we have silhouettes?
    if (is.null(x$sil))
      options$silhouettes <- FALSE

    if (options$reverse_columns)
      m <- m[, ncol(m):1]

    if (!options$silhouettes) {
      pushViewport(viewport(
        layout = grid.layout(
          6,
          3,
          widths = unit.c(
            unit(2, "lines"),
            # space
            unit(1, "snpc") - unit(7, "lines"),
            # image
            unit(2, "lines")                        # space
          ),
          heights = unit.c(
            unit(2, "lines"),
            # title
            unit(1, "lines"),
            # space
            unit(1, "snpc") - unit(7, "lines"),
            # image
            unit(1, "lines"),
            # space
            unit(1, "lines"),
            # colorkey
            unit(2, "lines")                        # space
          )
        ),
        gp = options$gp
      ))

      main_vp     <- viewport(
        layout.pos.col = 2,
        layout.pos.row = 1,
        name = "main"
      )
      image_vp    <-
        viewport(layout.pos.col = 2, layout.pos.row = 3)
      colorkey_vp <- viewport(
        layout.pos.col = 2,
        layout.pos.row = 5,
        name = "colorkey"
      )

    } else{
      ## with silhouettes
      pushViewport(viewport(
        layout = grid.layout(
          6,
          5,
          widths = unit.c(
            unit(2, "lines"),
            # space
            unit(0.7, "snpc") - unit(2.5, "lines"),
            # image
            unit(1, "lines"),
            # space
            unit(0.3, "snpc") - unit(2.5, "lines"),
            # sil
            unit(2, "lines")                        # space
          ),
          heights = unit.c(
            unit(2, "lines"),
            # title
            unit(2, "lines"),
            # space
            unit(0.7, "snpc") - unit(2.5, "lines"),
            # image
            unit(1, "lines"),
            # space
            unit(1, "lines"),
            # colorkey
            unit(2, "lines")                        # space
          )
        ),
        gp = options$gp
      ))

      main_vp     <-
        viewport(
          layout.pos.col = 2:4,
          layout.pos.row = 1,
          name = "main"
        )
      image_vp    <-
        viewport(layout.pos.col = 2,   layout.pos.row = 3)
      sil_vp  <- viewport(
        layout.pos.col = 4,
        layout.pos.row = 3,
        name = "sil"
      )
      colorkey_vp <-
        viewport(
          layout.pos.col = 2,
          layout.pos.row = 5,
          name = "colorkey"
        )

    }

    ## main
    pushViewport(main_vp)
    grid.text(options$main, gp = gpar(cex = 1.3, fontface = "bold"))
    upViewport(1)

    ## silhouette
    if (options$silhouettes) {
      ## get and reorder silhouettes
      s <- x$sil[, "sil_width"]

      pushViewport(sil_vp)
      .grid_barplot_horiz(s,
        xlab = "Silhouette width",
        gp_bars = gpar(fill = "lightgrey", col = 0))
      upViewport(1)

    }

    ## image
    if (is.null(options$zlim))
      options$zlim <- range(m, na.rm = TRUE)
    if (!is.null(options$threshold))
      m[m > options$threshold] <- NA

    pushViewport(image_vp)
    .grid_image(m, col = options$col, zlim = options$zlim)

    ## add labels?
    if (options$axes == "auto" && nrow(m) > 25)
      options$axes <- "none"
    if (options$axes != "none") {
      downViewport("image")
      #grid.text(colnames(m), y = unit(-1, "lines"),
      #  x=unit(1:ncol(m), "native"), rot=90, just="right")

      grid.text(
        rownames(m),
        x = unit(1, "npc") + unit(1, "lines"),
        y = unit(1:nrow(m), "native"),
        just = "left",
        gp = options$gp_labels
      )
      upViewport(1)
    }

    upViewport(1)

    ## color key?
    if (options$key) {
      pushViewport(colorkey_vp)
      .grid_colorkey(options$zlim,
        col = options$col,
        threshold = options$threshold)
      upViewport(1)
    }

    ## plot cluster borders if we have labels and order
    if (!is.null(labels)) {
      labels_unique_y		<- labels_unique
      cluster_width_y		<- (tabulate(labels)[labels_unique])
      cluster_cuts_y		<- cumsum(cluster_width_y)
      cluster_center_y	<- cluster_cuts_y - cluster_width_y / 2

      if (options$reverse_columns) {
        labels_unique_x	<- rev(labels_unique)
        cluster_width_x <- (tabulate(labels)[labels_unique_x])
        cluster_cuts_x  <- cumsum(cluster_width_x)
        cluster_center_x <- cluster_cuts_x - cluster_width_x / 2
      } else{
        labels_unique_x <- labels_unique_y
        cluster_width_x <- cluster_width_y
        cluster_cuts_x <- cluster_cuts_y
        cluster_center_x <- cluster_center_y
      }

      if (options$cluster_labels) {
        seekViewport("image")

        ## above the plot
        grid.text(
          labels_unique_x,
          x = cluster_center_x,
          y = unit(1, "npc") + unit(1, "lines"),
          default.units = "native",
          gp = options$gp_labels
        )
        ## left of the plot
        grid.text(
          labels_unique_y,
          x = unit(-1, "lines"),
          y = cluster_center_y,
          default.units = "native",
          gp = options$gp_labels
        )
        upViewport(2)
      }

      if (options$cluster_lines) {
        ## draw lines separating the clusters
        #cluster_cuts <- cluster_cuts[-length(cluster_cuts)]
        ## remove last line

        seekViewport("image")
        for (i in 1:(k - 1)) {
          grid.lines(
            #x = c(0, dim),
            x = c(0.5, dim + 0.5),
            y = cluster_cuts_y[i] + .5,
            default.units = "native",
            gp = options$gp_lines
          )

          grid.lines(
            x = cluster_cuts_x[i] + .5,
            #y = c(0, dim),
            y = c(0.5, dim + 0.5),
            default.units = "native",
            gp = options$gp_lines
          )

        }
        upViewport(2)

      }
    }

    if (options$pop)
      popViewport(1)
    else
      upViewport(1)
  }


## print for reordered_cluster_dissimilarity_matrix
#' @rdname dissplot
#' @export
print.reordered_cluster_dissimilarity_matrix <-
  function(x, ...)
  {
    d <- attr(x$x_reordered, "Size")
    k <- if (!is.null(x$k))
      x$k
    else
      NA

    cat(gettextf("object of class '%s'\n", class(x)))
    cat("matrix dimensions:", d, "x", d, "\n")
    cat(gettextf("dissimilarity measure: '%s'\n", x$diss_measure))
    cat("number of clusters k:", k, "\n")
    if (!is.null(x$k)) {
      cat("\ncluster description\n")
      print(x$description)
    }

    cat("\n")
    cat("used seriation methods\n")
    cat(gettextf("inter-cluster: '%s'\n", x$seriation_methods$inter))
    cat(gettextf("intra-cluster: '%s'\n", x$seriation_methods$intra))

    cat("\n")
    cat(gettextf(
      "dissimilarity aggregation method: '%s'\n",
      x$aggregation_method
    ))

    invisible(x)
  }

## inter and intra cluster dissimilarity matrix from
## a dissimilarity matrix plus labels
.cluster_dissimilarity <-
  function(x,
    labels,
    method = c("avg", "min", "max",
      "Hausdorff")) {
    method <- match.arg(method)
    ## FIXME: Implement Hausdorff

    linkage <- if (method == "avg")
      mean
    else if (method == "min")
      min
    else if (method == "max")
      max
    else if (method == "Hausdorff")
      .hausdorff
    else
      stop("Unknown method.")

    if (class(x) != "matrix")
      x <- as.matrix(x)


    ## kill self-dissimilarities (which are always 0)
    diag(x) <- NA

    k           <- length(unique(labels))
    diss_matrix <- matrix(nrow = k, ncol = k)

    ## calculate avg. dissimilarity between clusters
    for (i in 1:k) {
      slice <- x[labels == i, , drop = FALSE]
      for (j in 1:i) {
        block <- slice[, labels == j, drop = FALSE]

        val <- linkage(block, na.rm = TRUE)

        ## fix for clusters of size 1
        if (is.nan(val))
          val <- 0

        diss_matrix[i, j] <- val
        diss_matrix[j, i] <- val
      }
    }

    attr(diss_matrix, "method") <- method
    diss_matrix
  }

## implement Hausdorff distance between two sets from a dissimilarity matrix
##d_H = max{sup_x\inX inf_y\inY d(x,y), sup_y\inY inf_x\inX d(x,y)}
.hausdorff <- function(block, na.rm = TRUE)
  max(apply(block, MARGIN = 1, min, na.rm = na.rm),
    apply(block, MARGIN = 2, min, na.rm = na.rm))
