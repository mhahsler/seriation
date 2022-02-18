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

#' Visual Analysis for Cluster Tendency Assessment (VAT/iVAT)
#'
#' Implements Visual Analysis for Cluster Tendency Assessment (VAT; Bezdek and
#' Hathaway, 2002) and Improved Visual Analysis for Cluster Tendency Assessment
#' (iVAT; Wang et al, 2010).
#'
#' `path_dist()` redefines the distance between two objects as the minimum
#' over the largest distances in all possible paths between the objects as used
#' for iVAT.
#'
#' @family plots
#'
#' @param x a \code{dist} object.
#' @param upper_tri,lower_tri a logical indicating whether to show the upper or
#' lower triangle of the VAT matrix.
#' @param ... further arguments are passed on to \code{\link{pimage}} for the
#' regular plots and \code{\link{ggpimage}} for the ggplot2 plots.
#' @returns Nothing.
#'
#' @author Michael Hahsler
#' @references Bezdek, J.C. and Hathaway, R.J. (2002): VAT: a tool for visual
#' assessment of (cluster) tendency. \emph{Proceedings of the 2002
#' International Joint Conference on Neural Networks (IJCNN '02)}, Volume: 3,
#' 2225--2230.
#'
#' Havens, T.C. and Bezdek, J.C. (2012): An Efficient Formulation of the
#' Improved Visual Assessment of Cluster Tendency (iVAT) Algorithm, \emph{IEEE
#' Transactions on Knowledge and Data Engineering,} \bold{24}(5), 813--822.
#'
#' Wang L., U.T.V. Nguyen, J.C. Bezdek, C.A. Leckie and K. Ramamohanarao
#' (2010): iVAT and aVAT: Enhanced Visual Analysis for Cluster Tendency
#' Assessment, \emph{Proceedings of the PAKDD 2010, Part I, LNAI 6118,} 16--27.
#' @keywords cluster manip
#' @examples
#' ## lines data set from Havens and Bezdek (2011)
#' x <- create_lines_data(250)
#' plot(x, xlim=c(-5,5), ylim=c(-3,3), cex=.2)
#' d <- dist(x)
#'
#' ## create regular VAT
#' VAT(d, main = "VAT for Lines")
#' ## same as: pimage(d, seriate(d, "VAT"))
#'
#' ## ggplot2 version
#' if (require("ggplot2")) {
#'   ggVAT(d) + labs(title = "VAT")
#' }
#'
#' ## create iVAT which shows visually the three lines
#' iVAT(d, main = "iVAT for Lines")
#' ## same as:
#' ## d_path <- path_dist(d)
#' ## pimage(d_path, seriate(d_path, "VAT for Lines"))
#'
#' ## ggplot2 version
#' if (require("ggplot2")) {
#'   ggiVAT(d) + labs(title = "iVAT for Lines")
#' }
#'
#' ## compare with dissplot (shows banded structures and relationship between
#' ## center line and the two outer lines)
#' dissplot(d, method = "OLO_single", main = "Dissplot for Lines", col = bluered(100, bias = .5))
#'
#' ## compare with optimally reordered heatmap
#' hmap(d, method = "OLO_single", main = "Heatmap for Lines (opt. leaf ordering)",
#'   col = bluered(100, bias = .5))
#' @export
VAT <- function(x,
  upper_tri = TRUE,
  lower_tri = TRUE,
  ...) {
  if (!inherits(x, "dist"))
    stop("x needs to be of class 'dist'!")
  pimage(x,
    seriate(x, "VAT"),
    upper_tri = upper_tri,
    lower_tri = lower_tri,
    ...)
}

#' @rdname VAT
#' @export
iVAT <- function(x,
  upper_tri = TRUE,
  lower_tri = TRUE,
  ...) {
  if (!inherits(x, "dist"))
    stop("x needs to be of class 'dist'!")
  x <- path_dist(x)
  pimage(x,
    seriate(x, "VAT"),
    upper_tri = upper_tri,
    lower_tri = lower_tri,
    ...)
}


## calculate path distance from iVAT using a modified version fo Floyd's alg.
## d_ij = smallest value of the largest values of all possible paths between i and j

#' @rdname VAT
#' @export
path_dist <- function(x) {
  #A <- as.matrix(x)
  #n <- nrow(A)
  #for(k in 1:n)
  # for(i in 1:n)
  #   for(j in 1:n)
  #      if(max(A[i,k], A[k,j]) < A[i,j]) A[i,j] <- max(A[i,k], A[k,j])
  #d <- as.dist(A)

  ## make C call
  m <- as.matrix(x)

  if (any(is.na(m)))
    stop("NAs not allowed in x.")
  if (any(m < 0))
    stop("Negative values not allowed in x.")
  mode(m) <- "double"

  ## replace Inf with large number
  m[is.infinite(m)] <- .Machine$double.xmax

  if (any(m < 0))
    stop("Negative values not allowed in x.")

  m <- .Call("pathdist_floyd", m, PACKAGE = "seriation")
  as.dist(m)
}
