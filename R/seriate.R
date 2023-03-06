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

#' Seriate Dissimilarity Matrices, Matrices or Arrays
#'
#' Tries to find an linear order for objects using data in form of a
#' dissimilarity matrix (two-way one mode data), a data matrix (two-way
#' two-mode data) or a data array (k-way k-mode data). The order can then be
#' used to reorder the dissimilarity matrix/data matrix using
#' [permute()].
#'
#'
#' Seriation methods are managed via a registry. See
#' [list_seriation_methods()] for help. In the following, we discuss only the
#' built-in methods that are registered automatically by the package \pkg{seriation}.
#'
#' Many seriation methods (heuristically) optimize (minimize or maximize) an
#' objective function. The value of the function for a given seriation can be
#' calculated using [criterion()]. In this manual page, we
#' include the measure which is optimized by each method using **bold font**.
#' If no measure is mentioned, then the measure does not directly optimize a measure.
#' A definition of the measures can be found in the [criterion()] manual page.
#'
#' **Seriation methods for distance matrices (dist)**
#'
#' One-mode two-way data has to be provided as a dist object (not
#' as a symmetric matrix). Similarities have to be transformed into
#' dissimilarities. Currently, the following methods are implemented (for a
#' more detailed description and an experimental comparison see
#' [Hahsler (2017)](https://michael.hahsler.net/research/paper/EJOR_seriation_2016.pdf):
#'
#' - "ARSA" Anti-Robinson seriation by simulated
#'   annealing to minimize the **linear seriation criterion** (simulated
#'   annealing initialization used in Brusco et al 2008).
#'
#' - "BBURCG" Anti-Robinson seriation by branch-and-bound to
#'   minimize the **unweighted gradient measure** (Brusco and Stahl 2005).
#'   This is only feasible for a relatively small number of objects.
#'
#' - "BBWRCG" Anti-Robinson seriation by branch-and-bound to
#'   minimize the **weighted gradient measure** (Brusco and Stahl 2005). This
#'   is only feasible for a relatively small number of objects.
#'
#' - "TSP" Traveling salesperson problem solver to minimize the
#'   **Hamiltonian path length**. The solvers in \pkg{TSP} are used (see
#'   [TSP::solve_TSP()]). The solver method can be passed on via the `control`
#'   argument, e.g. `control = list(method = "two_opt")`. Default is the est
#'   of 10 runs of arbitrary insertion heuristic with 2-opt improvement.
#'
#'   Since a tour returned by a TSP solver is a connected circle and we are
#'   looking for a path representing a linear order, we need to find the best
#'   cutting point.  Climer and Zhang (2006) suggest to add a dummy city with
#'   equal distance to each other city before generating the tour. The place of
#'   this dummy city in an optimal tour with minimal length is the best cutting
#'   point (it lies between the most distant cities).
#'
#' - "R2E" Rank-two ellipse seriation (Chen 2002).
#'
#'   This method starts with generating a sequence of correlation matrices
#'   \eqn{R^1, R^2, \ldots}. \eqn{R^1} is the correlation matrix of the original
#'   distance matrix \eqn{D} (supplied to the function as `x`), and
#'   \deqn{R^{n+1} = \phi R^n,} where \eqn{\phi} calculates the correlation
#'   matrix.
#'
#'   The rank of the matrix \eqn{R^n} falls with increasing \eqn{n}. The first
#'   \eqn{R^n} in the sequence which has a rank of 2 is found. Projecting all
#'   points in this matrix on the first two eigenvectors, all points fall on an
#'   ellipse. The order of the points on this ellipse is the resulting order.
#'
#'   The ellipse can be cut at the two interception points (top or bottom) of the
#'   vertical axis with the ellipse. In this implementation the top most cutting
#'   point is used.
#'
#' - "MDS", "MDS_metric", "MDS_nonmetric", "MDS_angle" Multidimensional scaling (MDS).
#'
#'   Use multidimensional scaling techniques to find an linear order by
#'   minimizing **stress**. Note MDS algorithms used for a single dimension
#'   tend to end up in local optima and unidimensional scaling (see Maier and De
#'   Leeuw, 2015) would be more appropriate. However, generally, ordering along
#'   the first component of MDS provides good results.
#'
#'   `control` parameters:
#'     - `method`: One of `"cmdscale"`, `"isoMDS"` or `"sammon"`. `"cmdscale"` performs metric
#'        MDS using [stats::cmdscale()]. Non-metric MDS methods `"isoMDS"` and `"sammon"`
#'        are preformed using [MASS::isoMDS()].
#'
#'   By default, metric MDS is used ([stats::cmdscale()]). In case of
#'   of general dissimilarities, non-metric MDS can be used. The method can be specified
#'   as the element `method` (`"cmdscale"`, `"isoMDS"` or `"sammon"`) in `control`.
#'
#'   For convenience, seriation methods `"MDS_metric"` performs [cmdscale()] and
#'   `"MDS_nonmetric"` performs [MASS::isoMDS()].
#'
#'   `"MDS_angle"` projects the data on the first two components found by
#'   MDS and then orders by the angle in this space. The order is split by the
#'   larges gap between adjacent angles. A similar method was used for ordering
#'   correlation matrices by Friendly (2002).
#'
#' - "HC", "HC_single", "HC_complete", "HC_average", "HC_ward" Hierarchical clustering.
#'
#'   Using the order of the leaf nodes in a dendrogram obtained by hierarchical
#'   clustering can be used as a very simple seriation technique. This method
#'   applies hierarchical clustering ([hclust()]) to `x`. The clustering
#'   method can be given using a `"method"` element in the `control`
#'   list. If omitted, the default `"average"` is used.
#'
#'   For convenience the other methods are provided as shortcuts.
#'
#' - "GW" Hierarchical clustering (Gruvaeus and Wainer, 1972).
#'
#'   The methods start with a dendrogram created by [hclust()]. As the
#'   `"method"` element in the `control` list a clustering method
#'   (default `"average"`) can be specified. Alternatively, an [hclust]
#'   object can be supplied using an element named `"hclust"`.
#'
#'   A dendrogram (binary tree) has \eqn{2^{n-1}} internal nodes (subtrees) and
#'   the same number of leaf orderings. That is, at each internal node the left
#'   and right subtree (or leaves) can be swapped, or, in terms of a dendrogram,
#'   be flipped. The leaf-node reordering to minimize
#'   **Hamiltonian path length (restricted)**.
#'
#'   Method `"GW"` uses an algorithm developed by Gruvaeus and Wainer (1972)
#'   as implemented [gclus::reorder.hclust()] (Hurley 2004).  The clusters are
#'   ordered at each level so that the objects at the edge of each cluster are
#'   adjacent to that object outside the cluster to which it is nearest. The
#'   method produces an unique order.
#'
#'   For convenience `"GW_single"`, `"GW_average"`,
#'   `"GW_complete"`, and `"GW_ward"` are provided.
#'
#' - "OLO" Optimal leaf ordering (Bar-Joseph et al., 2001).
#'
#'   Also starts with a dendrogram and
#'   produces an optimal leaf ordering with respect to the minimizing the sum of
#'   the distances along the (Hamiltonian) path connecting the leaves in the
#'   given order. The time complexity of the algorithm is \eqn{O(n^3)}. Note that
#'   non-finite distance values are not allowed.
#'
#'   For convenience `"OLO_single"`,
#'   `"OLO_average"`, `"OLO_complete"`, and `"OLO_ward"` are provided.
#'
#' - "VAT" Visual Assessment of (Clustering) Tendency (Bezdek and Hathaway (2002)).
#'
#'   Creates an order based on Prim's algorithm for finding a minimum spanning
#'   tree (MST) in a weighted connected graph representing the distance matrix.
#'   The order is given by the order in which the nodes (objects) are added to
#'   the MST.
#'
#' - "SA" Simulated Annealing for diverse criterion measures.
#'
#'   Implement simulated annealing similar to the ARSA method, however, it works
#'   for any criterion measure defined in \pkg{seriation}. By default the
#'   algorithm optimizes for raw gradient measure and is warm started with the
#'   result of spectral seriation (2-Sum problem) since Hahsler (2017) shows that
#'   2-Sum solutions are similar to solutions for the gradient measure.
#'
#'   Several popular local neighborhood functions are
#'   provided an new can be defined (see [LS]).
#'
#'   Note that this is an R implementation repeatedly calling criterion, and
#'   therefore is relatively slow.
#'
#' - "Spectral", "Spectral_norm" Spectral seriation (Ding and He 2004).
#'
#'   Spectral seriation uses a relaxation to minimize the **2-Sum Problem**
#'   (Barnard, Pothen, and Simon, 1993). It uses the order of the Fiedler vector
#'   of the similarity matrix's (normalized) Laplacian.
#'
#'   Spectral seriation gives a good trade-off between seriation quality, speed
#'   and scalability (see Hahsler, 2017).
#'
#' - "SPIN_STS", "SPIN_NH" Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005).
#'
#'   Given a weight matrix \eqn{W}, the algorithms try to
#'   minimize the energy for a permutation (matrix \eqn{P}) given by \deqn{F(P) =
#'   tr(PDP^TW),} where \eqn{tr} denotes the matrix trace.
#'
#'   `"SPIN_STS"` implements the Side-to-Side algorithm which tries to push
#'   out large distance values. The default weight matrix suggested in the paper
#'   with \eqn{W=XX^T} and \eqn{X_i=i-(n+1)/2} is used. We run the algorithm from
#'   `step` (25) iteration and restart the algorithm `nstart` (10) with
#'   random initial permutations (default values in parentheses). Via
#'   `control` the parameters `step`, `nstart`, `X` and
#'   `verbose`.
#'
#'   `"SPIN_NH"` implements the neighborhood algorithm (concentrate low
#'   distance values around the diagonal) with a Gaussian weight matrix
#'   \eqn{W_{ij} = exp(-(i-j)^2/n\sigma)}, where \eqn{n} is the size of the
#'   dissimilarity matrix and \eqn{\sigma} is the variance around the diagonal
#'   that control the influence of global (large \eqn{\sigma}) or local (small
#'   \eqn{\sigma}) structure.
#'
#'   We use the heuristic suggested in the paper for the linear assignment
#'   problem. We do not terminate as indicated in the algorithm, but run all the
#'   iterations since the heuristic does not guarantee that the energy is
#'   strictly decreasing. We also implement the heuristic "annealing" scheme
#'   where \eqn{\sigma} is successively reduced. The parameters in `control`
#'   are `sigma` which can be a single value or a decreasing sequence
#'   (default: 20 to 1 in 10 steps) and `step` which defines how many update
#'   steps are performed before for each value of `alpha`. Via
#'   `W_function` a custom function to create \eqn{W} with the function
#'   signature `function(n, sigma, verbose)` can be specified. The parameter
#'   `verbose` can be used to display progress information.
#'
#' - "QAP_LS", "QAP_2SUM", "QAP_BAR", "QAP_Inertia" Quadratic assignment problem
#'   formulations for seriation using a simulated annealing solver.
#'
#'   These methods minimize the
#'   **Linear Seriation Problem** (LS) formulation (Hubert and Schultz 1976),
#'   the **2-Sum Problem** formulation (Barnard, Pothen, and Simon 1993), the
#'   **banded anti-Robinson form** (BAR) or the **inertia criterion.**
#'
#'   `control` parameters are passed on to [qap::qap()].
#'   An important parameter is `rep` to return the best result out of the
#'   given number of repetitions with random restarts. Default is 1, but bigger
#'   numbers result in better and more stable results.
#'
#' - "GA" Use a genetic algorithm to optimize for various criteria.
#'
#'   The GA code has to be first registered. A detailed description can
#'   be found in the manual page for [register_GA()].
#'
#' - "DendSer" Use heuristic dendrogram seriation to optimize for various criteria.
#'
#'   The DendSer code has to be first registered. A
#'   detailed description can be found in the manual page for
#'   [register_DendSer()].
#'
#' - "Identity" Produces an identity permutation.
#'
#' - "Random"  Produces a random permutation.
#'
#'
#' **Seriation methods for matrices (matrix or data.frame)**
#'
#' Two-mode two-way data are general matrices.
#' Some methods also require that the matrix is positive. Data frames are just a
#' different representation of a matrix and all seriation methods for matrix can
#' be also used for data frames. The default method for data frames is heatmap
#' seriation which calculates distances between rows and between columns and
#' then applies seriation on these using hierarchical clustering and optimal
#' leaf ordering (OLO).
#'
#' Currently the
#' following methods are implemented for matrix:
#'
#' - "BEA" Bond Energy Algorithm (BEA; McCormick 1972).
#'
#'   The algorithm tries to maximize the **Measure of Effectiveness.** of a
#'   non-negative matrix. Due to the definition of this measure, the tasks of
#'   ordering rows and columns is separable and can be solved independently.
#'
#'   A row is arbitrarily placed; then rows are positioned one by one. When this
#'   is completed, the columns are treated similarly. The overall procedure
#'   amounts to two approximate traveling salesperson problems (TSP), one on the
#'   rows and one on the columns. The so-called `best insertion strategy' is
#'   used: rows (or columns) are inserted into the current permuted list of rows
#'   (or columns). Several consecutive runs of the algorithm might improve the
#'   energy.
#'
#'   Note that Arabie and Hubert (1990) question its use with non-binary data if
#'   the objective is to find a seriation or one-dimensional orderings of rows
#'   and columns.
#'
#'   The BEA code used in this package was implemented by Fionn Murtagh.
#'
#'   `control` parameter:
#'     - `"rep"`: the number of runs can be specified.
#'        The results of the best run will be returned.
#'
#' - "BEA_TSP" Use a TSP to optimize the **Measure of Effectiveness** (Lenstra 1974).
#'
#'   `control` parameter:
#'      - `"method"`: a TSP solver method (see [TSP::solve_TSP()]).
#'
#' - "CA" Correspondence analysis for a table/matrix of frequencies.
#'
#'   This function is designed to help simplify a mosaic plot or other displays of a
#'   matrix of frequencies.  It calculates a correspondence analysis of the matrix and
#'   an order for rows and columns according to the scores on a correspondence analysis dimension.
#'
#'   `control` parameters:
#'     - `"dim"`: CA dimension used for reordering.
#'     - `"ca_param"`: List with parameters for the call to [ca::ca()].
#'
#' - "Heatmap" Heatmap seriation
#'
#'   Calculates distances between
#'   rows and between columns and then applies seriation on these using
#'   hierarchical clustering and optimal leaf ordering (method `"OLO"` for distance matrices).
#'
#' - "PCA" Order by the first principal component.
#'
#'   Uses the projection of the data on its first principal component to
#'   determine the order.
#'
#'   Note that for a distance matrix calculated from `x` with Euclidean
#'   distance, this methods minimizes the least square criterion.
#'
#' - "PCA_angle" Order using the first two principal components.
#'
#'   Projects the data on the first two principal components
#'   and then orders by the angle in this space. The order is split by the larges
#'   gap between adjacent angles. A similar method was used for ordering
#'   correlation matrices by Friendly (2002).
#'
#' - "Identity" Produces an identity permutation.
#'
#' - "Random" Produces a random permutation.
#'
#' For **general arrays** no built-in methods are currently available.
#'
#' @family seriation
#'
#' @param x the data.
#' @param method a character string with the name of the seriation method
#' (default: varies by data type).
#' @param control a list of control options passed on to the seriation
#' algorithm.
#' @param margin a vector giving the margin indices (dimensions) to be
#' seriated. For example, for a matrix, `1` indicates rows, `2`
#' indicates columns, `c(1,2)` indicates rows and columns. Unseriated margins return
#' a identity seriation order.
#' @param ... further arguments are added to the `control` list.
#'
#' @return Returns an object of class [ser_permutation].
#'
#' @author Michael Hahsler
#'
#' @references Arabie, P. and L.J. Hubert (1990): The bond energy algorithm
#' revisited, \emph{IEEE Transactions on Systems, Man, and Cybernetics,}
#' \bold{20}(1), 268--274.
#' \doi{10.1109/21.47829}
#'
#' Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola. (2001): Fast
#' Optimal Leaf Ordering for Hierarchical Clustering. \emph{Bioinformatics,}
#' \bold{17}(1), 22--29.
#' \doi{10.1093/bioinformatics/17.suppl_1.S22}
#'
#' Barnard, S. T., A. Pothen, and H. D. Simon (1993): A Spectral Algorithm for
#' Envelope Reduction of Sparse Matrices. \emph{In Proceedings of the 1993
#' ACM/IEEE Conference on Supercomputing,} 493--502. Supercomputing '93. New
#' York, NY, USA: ACM. \url{https://ieeexplore.ieee.org/document/1263497}
#'
#' Bezdek, J.C. and Hathaway, R.J. (2002): VAT: a tool for visual assessment of
#' (cluster) tendency. \emph{Proceedings of the 2002 International Joint
#' Conference on Neural Networks (IJCNN '02)}, Volume: 3, 2225--2230.
#' \doi{10.1109/IJCNN.2002.1007487}
#'
#' Brusco, M., Koehn, H.F., and Stahl, S. (2008): Heuristic Implementation of
#' Dynamic Programming for Matrix Permutation Problems in Combinatorial Data
#' Analysis. \emph{Psychometrika,} \bold{73}(3), 503--522.
#' \doi{10.1007/s11336-007-9049-5}
#'
#' Brusco, M., and Stahl, S. (2005): \emph{Branch-and-Bound Applications in
#' Combinatorial Data Analysis.} New York: Springer.
#' \doi{10.1007/0-387-28810-4}
#'
#' Chen, C. H. (2002): Generalized Association Plots: Information Visualization
#' via Iteratively Generated Correlation Matrices. \emph{Statistica Sinica,}
#' \bold{12}(1), 7--29.
#'
#' Ding, C. and Xiaofeng He (2004): Linearized cluster assignment via spectral
#' ordering. \emph{Proceedings of the Twenty-first International Conference on
#' Machine learning (ICML '04)}.
#' \doi{10.1145/1015330.1015407}
#'
#' Climer, S. and Xiongnu Zhang (2006): Rearrangement Clustering: Pitfalls,
#' Remedies, and Applications, \emph{Journal of Machine Learning Research,}
#' \bold{7}(Jun), 919--943.
#'
#' Friendly, M. (2002): Corrgrams: Exploratory Displays for Correlation
#' Matrices. \emph{The American Statistician}, \bold{56}(4), 316--324.
#' \doi{10.1198/000313002533}
#'
#' Gruvaeus, G. and Wainer, H. (1972): Two Additions to Hierarchical Cluster
#' Analysis, \emph{British Journal of Mathematical and Statistical Psychology,}
#' \bold{25}, 200--206.
#' \doi{10.1111/j.2044-8317.1972.tb00491.x}
#'
#' Hahsler, M. (2017): An experimental comparison of seriation methods for
#' one-mode two-way data. \emph{European Journal of Operational Research,}
#' \bold{257}, 133--143.
#' \doi{10.1016/j.ejor.2016.08.066}
#'
#' Hubert, Lawrence, and James Schultz (1976): Quadratic Assignment as a
#' General Data Analysis Strategy. \emph{British Journal of Mathematical and
#' Statistical Psychology} \bold{29}(2). Blackwell Publishing Ltd. 190--241.
#' \doi{10.1111/j.2044-8317.1976.tb00714.x}
#'
#' Hurley, Catherine B. (2004): Clustering Visualizations of Multidimensional
#' Data. \emph{Journal of Computational and Graphical Statistics,}
#' \bold{13}(4), 788--806.
#' \doi{10.1198/106186004X12425}
#'
#' Lenstra, J.K (1974): Clustering a Data Array and the Traveling-Salesman
#' Problem, \emph{Operations Research,} \bold{22}(2) 413--414.
#' \doi{10.1287/opre.22.2.413}
#'
#' Mair P., De Leeuw J. (2015). Unidimensional scaling. In \emph{Wiley
#' StatsRef: Statistics Reference Online,} Wiley, New York.
#' \doi{10.1002/9781118445112.stat06462.pub2}
#'
#' McCormick, W.T., P.J. Schweitzer and T.W. White (1972): Problem
#' decomposition and data reorganization by a clustering technique,
#' \emph{Operations Research,} \bold{20}(5), 993--1009.
#' \doi{10.1287/opre.20.5.993}
#'
#' Tsafrir, D., Tsafrir, I., Ein-Dor, L., Zuk, O., Notterman, D.A. and Domany,
#' E. (2005): Sorting points into neighborhoods (SPIN): data analysis and
#' visualization by ordering distance matrices, \emph{Bioinformatics,}
#' \bold{21}(10) 2301--8.
#' \doi{10.1093/bioinformatics/bti329}
#' @keywords optimize cluster
#' @examples
#' # Show available seriation methods (for dist and matrix)
#' list_seriation_methods()
#'
#' ### Seriate as distance matrix (for the iris dataset)
#' data("iris")
#' x <- as.matrix(iris[-5])
#' x <- x[sample(1:nrow(x)), ]
#' d <- dist(x)
#'
#' order <- seriate(d)
#' order
#'
#' pimage(d, main = "Distances (Random Order)")
#' pimage(d, order, main = "Distances (Reordered)")
#'
#' # Compare seriation quality
#' rbind(
#'         random = criterion(d),
#'         reordered = criterion(d, order)
#'      )
#'
#' # Reorder the distance matrix
#' d_reordered <-  permute(d, order)
#' pimage(d_reordered, main = "Distances (Reordered)")
#'
#'
#' ### Seriate a matrix
#' data("iris")
#' x <- as.matrix(iris[-5])
#'
#' # To make the variables comparable, we scale the data
#' x <- scale(x, center = FALSE)
#'
#' # The iris flowers are ordered by species in the data set
#' pimage(x, main = "original data", prop = FALSE)
#' criterion(x)
#'
#' # Apply some methods
#' order <- seriate(x, method = "BEA_TSP")
#' pimage(x, order, main = "TSP to optimize ME", prop = FALSE)
#' criterion(x, order)
#'
#' order <- seriate(x, method = "PCA")
#' pimage(x, order, main = "First principal component", prop = FALSE)
#' criterion(x, order)
#'
#' order <- seriate(x, method = "heatmap")
#' pimage(x, order, main = "Heatmap seriation", prop = FALSE)
#' criterion(x, order)
#'
#' # reorder the matrix
#' x_reordered <- permute(x, order)
#'
#' # create a heatmap seriation manually by calculating
#' # distances between rows and between columns
#' order <- c(
#'     seriate(dist(x), method = "OLO"),
#'     seriate(dist(t(x)), method = "OLO")
#' )
#' pimage(x, order, main = "Heatmap seriation", prop = FALSE)
#' criterion(x, order)
#'
#' ### Seriate a correlation matrix
#' corr <- cor(x)
#' pimage(corr, upper_tri = FALSE, main = "Correlation matrix")
#'
#' # we need to define a distance (we used d = sqrt(2(1 - r))) and
#' # then reorder the matrix (rows and columns).
#' d <- as.dist(sqrt(2 * (1 - corr)))
#' o <- seriate(d)
#' corr_reordered <- permute(corr, order = c(o, o))
#' pimage(corr_reordered, upper_tri = FALSE, main = "Correlation matrix (reordered)")
#' @export
seriate <- function(x, ...)
  UseMethod("seriate")

#' @export
seriate.default <- function(x, ...)
  stop(gettextf("seriate not implemented for class '%s'.",
    class(x)))
