#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2017 Michael Hahsler, Christian Buchta and Kurt Hornik
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

#' Unidimensional Scaling from Seriation Results
#'
#' Calculates an (approximate) unidimensional scaling configuration given an order.
#'
#' This implementation uses the method describes in Maier and De Leeuw (2015) to calculate the
#' minimum stress configuration for a given (seriation) order by performing a 1D MDS fit.
#' If the 1D MDS fit does not preserve the given order perfectly, then a warning is
#' produced indicating
#' for how many positions order could not be preserved.
#' The seriation method which is consistent to uniscale is `"MDS_smacof"`
#' which needs to be registered with [`register_smacof()`].
#'
#'
#' The code is similar to `uniscale()` in \pkg{smacof} (de Leeuw, 2090), but scales to larger
#' datasets since it only uses the permutation given by  `order`.
#'
#' `MDS_stress` calculates the normalized stress of a configuration given by a seriation order.
#' If the order does not contain a configuration, then a minimum-stress configuration if calculates
#' for the given order.
#'
#' All distances are first normalized to an average distance of close to 1 using
#' \eqn{d_{ij} \frac{\sqrt{n(n-1)/2}}{\sqrt{\sum_{i<j}{d_{ij}}^2}}}.
#'
#' Some seriation methods produce a MDS configuration (a 1D or 2D embedding). `get_config()`
#' retrieved the configuration attribute from the `ser_permutation_vector`. `NULL`
#' is returned if the seriation did not produce a configuration.
#'
#' `plot_config()` plots 1D and 2D configurations. `...` is passed on
#'   to [`plot.default`] and accepts `col`, `labels`, etc.
#'
#' @param d a dissimilarity matrix.
#' @param order a precomputed permutation (configuration) order.
#' @param accept_reorder logical; accept a configuration that does not preserve
#'  the requested order. If `FALSE`, the initial configuration stored in `order`
#'  or, an equally spaced configuration is returned.
#' @param warn logical; produce a warning if the 1D MDS fit does not preserve the
#'  given order.
#' @param \dots additional arguments are passed on to the seriation method.
#' @return A vector with the fitted configuration.
#' @author Michael Hahsler with code from Patrick Mair (from \pkg{smacof}).
#' @references Mair P., De Leeuw J. (2015). Unidimensional scaling. In
#' \emph{Wiley StatsRef: Statistics Reference Online,} Wiley, New York.
#' \doi{10.1002/9781118445112.stat06462.pub2}
#'
#' Jan de Leeuw, Patrick Mair (2009). Multidimensional Scaling Using Majorization:
#' SMACOF in R. Journal of Statistical Software, 31(3), 1-30.
#' \doi{10.18637/jss.v031.i03}
#' @seealso [register_smacof()]
#' @keywords optimize
#' @examples
#' data(SupremeCourt)
#' d <- as.dist(SupremeCourt)
#' d
#'
#' # embedding-based methods return "configuration" attribute
#' # plot_config visualizes the configuration
#' o <- seriate(d, method = "sammon")
#' get_order(o)
#' plot_config(o)
#'
#' # the configuration (Note: objects are in the original order in d)
#' get_config(o)
#'
#' # angle methods return a 2D configuration
#' o <- seriate(d, method = "MDS_angle")
#' get_order(o)
#' get_config(o)
#' plot_config(o, )
#'
#'
#' # calculate a configuration for a seriation method that does not
#' # create a configuration
#' o <- seriate(d, method = "ARSA")
#' get_order(o)
#' get_config(o)
#'
#' # find the minimum-stress configuration for the ARSA order
#' sc <- uniscale(d, o)
#' sc
#'
#' plot_config(sc)
#' @export
uniscale <-
  function(d,
           order,
           accept_reorder = FALSE,
           warn = TRUE,
           ...) {
    d <- as.dist(d)
    order <- ser_permutation(order)

    init_config <- get_config(order)

    # we cannot use more than 1D
    if (is.matrix(init_config))
      init_config <- NULL

    if (is.null(init_config))
      init_config <- get_rank(order)
    n <- length(init_config)

    if (attr(d, "Size") != n)
      stop("size of d and length of order do not agree!")

    # we do not use weights
    w <- 1 - diag(n)
    delta <- as.matrix(.normDiss(d))

    v <- as.matrix(solve((diag(rowSums(
      w
    )) - w) + (1 / n)) - (1 / n))
    s <- sign(outer(init_config, init_config, "-"))

    t <- as.vector(v %*% rowSums(delta * w * s))

    # does the configuration preserve the order in o?
    mismatches <- sum(order(init_config) != order(t))
    if (mismatches > 0 && warn) {
      warning("Configutation order does not preserve given order! Mismatches: ",
              mismatches,
              " of ",
              n, " - returning initial configuration instead.")
    }

    if (!accept_reorder && mismatches > 0)
      t <- init_config

    #cat("init:\n")
    #print(names(init_config))
    #cat("d:\n")
    #print(labels(d))

    names(t) <- labels(d)
    t
  }

# normalize the distances to roughly n*(n-1) / 2 so the average distance
# is close to 1
.normDiss <- function (diss)
  diss / sqrt(sum(diss ^ 2, na.rm = TRUE)) *
  sqrt(length(diss))

#' @rdname uniscale
#' @param refit logical; forces to refit a minimum-stress MDS configuration,
#'  even if `order` contains a configuration.
#' @export
MDS_stress <- function(d, order, refit = TRUE, warn = FALSE) {

  d <- as.dist(d)
  o <- ser_permutation(order)

  emb <- get_config(o)
  if(is.null(emb) || refit)
    emb <- uniscale(d, o, warn = warn)

  d_emb <- dist(emb)

  d_emb <- .normDiss(d_emb)
  d <- .normDiss(d)

  sqrt(sum((d - d_emb)^2) / sum(d_emb^2))
}

set_criterion_method(
  "dist",
  "MDS_stress",
  MDS_stress,
  "Normalized stress of a configuration given by a seriation order",
  FALSE
)

#' @rdname uniscale
#' @param dim The dimension if `x` is a `ser_permutation` object.
#' @export
get_config <- function(x, dim = 1L, ...) {
  if (inherits(x, "ser_permutation"))
    x <- x[[dim]]

  if (inherits(x, "ser_permutation_vector"))
    x <- attr(x, "configuration")

  if(is.null(x))
    return(NULL)

  if (!(is.numeric(x) && ((is.vector(x) || is.matrix(x)))))
    stop("Unable to get configuration. Supply a ser_permutation.")

  x
}


#' @rdname uniscale
#' @param x a scaling returned by `uniscale()` or a
#'   `ser_permutation` with a configuration attribute.
#' @param main main plot label
#' @param pch print character
#' @param labels add the object names to the plot
#' @param pos label position for 2D plot (see [text()]).
#' @param cex label expansion factor.
#' @export
plot_config <- function (x,
                        main,
                        pch = 19,
                        labels = TRUE,
                        pos = 1,
                        cex = 1,
                        ...) {
  if (missing(main))
    main <- "Configuration"

  o <- get_order(x)
  x <- get_config(x)


  if (is.null(x))
    stop(
      "Permutation vector has no configuration attribute. Use uniscale() first to calcualte a configuration"
    )

  # 2D
  if (is.matrix(x)) {
    graphics::plot(x, pch = pch, main = main, ...)

    if (labels)
      graphics::text(x = x,
                     labels = rownames(x),
                     pos = pos,
                     cex = cex)
    graphics::lines(x[get_order(o), , drop = FALSE], col = "grey")

  } else{
    # 1D
    x <- drop(x)
    n <- length(x)
    plot(
      x,
      rep(0, n),
      axes = FALSE,
      ann = FALSE,
      pch = pch,
      type = "o",
      ylim = c(-0.2, 0.8),
      ...
    )
    title(main)

    labs <- names(x)
    if (is.null(labs))
      labs <- 1:n

    if (labels)
      text(x,
           rep(0, n) + 0.05,
           labs,
           srt = 90,
           cex = cex,
           adj = c(0, 0.5))
  }
}
