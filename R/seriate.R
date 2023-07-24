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
#' Tries to find a linear order for objects using data in the form of a
#' dissimilarity matrix (two-way one-mode data), a data matrix (two-way
#' two-mode data), or a data array (k-way k-mode data). The order can then be
#' used to reorder the dissimilarity matrix/data matrix using
#' [permute()].
#'
#' Seriation methods are managed via a registry. See
#' [list_seriation_methods()] for help. In the following, we focus on
#' discussing the
#' built-in methods that are registered automatically by the package \pkg{seriation}.
#'
#' The available control options, default settings, and
#' a description for each algorithm
#' can be retrieved using `get_seriation_method(name = "<seriation method>")`.
#' Some control parameters are also described in more detail below.
#'
#' Some methods are very slow, and progress can be printed using the control
#' parameter `verbose = TRUE`.
#'
#' Many seriation methods (heuristically) optimize (minimize or maximize) an
#' objective function often called seriation criterion.
#' The value of the seriation criterion for a given order can be
#' calculated using [criterion()]. In this manual page, we
#' include the criterion, which is optimized by each method using **bold font**.
#' If no criterion is mentioned, then the method does not directly optimize a criterion.
#' A definition of the different seriation criteria can be found on the [criterion()] manual page.
#'
#' **Seriation methods for distance matrices (dist)**
#'
#' One-mode two-way data must be provided as a dist object (not
#' a symmetric matrix). Similarities have to be transformed into
#' dissimilarities.
#' Seriation algorithms fall into different groups based on the approach.
#' In the following, we describe the currently implemented methods.
#' A list with all methods and the available parameters is available
#' [here](https://mhahsler.github.io/seriation/seriation_methods.html).
#' [Hahsler (2017)](https://michael.hahsler.net/research/paper/EJOR_seriation_2016.pdf)
#' for a more detailed description and an experimental comparison of the most
#' popular methods.
#'
#'
#' **Dendrogram leaf order**
#'
#'  These methods create a dendrogram using hierarchical clustering and then derive
#'  the seriation order from the leaf order in the dendrogram. Leaf reordering
#'  may be applied.
#'
#'  - **Hierarchical clustering:** `"HC"`, `"HC_single"`, `"HC_complete"`,
#'      `"HC_average"`, `"HC_ward"`
#'
#'       Uses the order of the leaf nodes in a dendrogram obtained by hierarchical
#'       clustering as a simple seriation technique. This method
#'       applies hierarchical clustering ([hclust()]) to `x`. The clustering
#'       method can be given using a `"linkage"` element in the `control`
#'       list. If omitted, the default `"complete"` is used.
#'       For convenience, the other methods are provided as shortcuts.
#'
#' - **Reordered by the Gruvaeus and Wainer heuristic:** `"GW"`, `"GW_single"`, `"GW_average"`,
#'   `"GW_complete"`, `"GW_ward"`  (Gruvaeus and Wainer, 1972)
#'
#'   Method `"GW"` uses an algorithm developed by Gruvaeus and Wainer (1972)
#'   as implemented [gclus::reorder.hclust()] (Hurley 2004).  The clusters are
#'   ordered at each level so that the objects at the edge of each cluster are
#'   adjacent to the nearest object outside the cluster. The
#'   method produces a unique order.
#'
#'     The methods start with a dendrogram created by [hclust()]. As the
#'     `"linkage"` element in the `control` list, a clustering method
#'     (default `"average"`) can be specified. Alternatively, an [hclust]
#'     object can be supplied using an element named `"hclust"`.
#'
#'     A dendrogram (binary tree) has \eqn{2^{n-1}} internal nodes (subtrees) and
#'     the same number of leaf orderings. That is, at each internal node, the left
#'     and right subtree (or leaves) can be swapped or, in terms of a dendrogram,
#'     be flipped. The leaf-node reordering to minimize
#'
#'     Minimizes the **Hamiltonian path length (restricted by the dendrogram)**.
#'
#' - **Reordered by optimal leaf ordering:** `"OLO"`, `"OLO_single"`,
#'   `"OLO_average"`, `"OLO_complete"`, `"OLO_ward"`  (Bar-Joseph et al., 2001)
#'
#'   Starts with a dendrogram and
#'   produces an optimal leaf ordering that minimizes the sum of
#'   the distances along the (Hamiltonian) path connecting the leaves in the
#'   given order. The algorithm's time complexity is \eqn{O(n^3)}. Note that
#'   non-finite distance values are not allowed.
#'
#'   Minimizes the **Hamiltonian path length (restricted by the dendrogram)**.
#'
#' - **Dendrogram seriation:** `"DendSer"` (Earle and Hurley, 2015)
#'
#'    Use heuristic dendrogram seriation to optimize for various criteria.
#'    The DendSer code has to be first registered. A
#'    detailed description can be found on the manual page for
#'    [register_DendSer()].
#'
#' **Dimensionality reduction**
#'
#' Find a seriation order by reducing the dimensionality to 1 dimension. This is typically
#' done by minimizing a stress measure or the reconstruction error.
#' Note that dimensionality reduction to a single dimension is a very
#' difficult discrete optimization problem.
#' For example, MDS algorithms used for a single dimension
#' tend to end up in local optima (see Maier and De Leeuw, 2015).
#' However, generally, ordering along a single component of MDS provides good results
#' sufficient for applications like visualization.
#'
#'   - **Classical metric multidimensional scaling:** `"MDS"`
#'
#'     Orders along the 1D classical metric multidimensional scaling.
#'     `control` parameters are passed on to [stats::cmdscale()].

#'   - **Isometric feature mapping:** `"isomap"` (Tenenbaum, 2000)
#'
#'     Orders along the 1D isometric feature mapping.
#'     `control` parameters are passed on to [vegan::isomap()]
#'
#'   - **Kruskal's non-metric multidimensional scaling:** `"isoMDS"`, `"monoMDS"`,
#'    `"metaMDS"` (Kruskal, 1964)
#'
#'      Orders along the 1D Kruskal's non-metric multidimensional scaling.
#'      Package `vegan` implements an alternative implementation called `monoMDS`
#'      and a version that uses random restarts for stability called `metaMDS`.
#'      `control` parameters are passed on to [MASS::isoMDS()], [vegan::monoMDS()] or [vegan::metaMDS()].
#'
#'   - **Sammon's non-linear mapping:** `"Sammon_mapping"` (Sammon, 1969)
#'
#'       Orders along the 1D Sammon's non-linear mapping.
#'       `control` parameters are passed on to [MASS::sammon()].
#'
#'
#'   - **Angular order of the first two eigenvectors:** `"MDS_angle"`
#'
#'       Finds a 2D configuration using MDS ([cmdscale()])
#'       to approximate the eigenvectors of the covariance matrix in the
#'       original data matrix.
#'       Orders by the angle in this space and splits the order by the
#'       larges gap between adjacent angles. A similar method was used by
#'       Friendly (2002) to order variables in correlation matrices
#'       by angles of first two eigenvectors.
#'
#'   - **Smacof:** `"MDS_smacof"` (de Leeuw and Mair, 2009)
#'
#'      Perform seriation using stress majorization with several transformation functions.
#'      This method has to be registered first using [`register_smacof()`].
#'
#' **Optimization**
#'
#' These methods try to optimize a seriation criterion directly, typically using a
#' heuristic approach.
#'
#' - **Anti-Robinson seriation by simulated annealing:** `"ARSA"` (Brusco et al 2008)
#'
#'   The algorithm automatically finds a suitable start temperature and calculates
#'   the needed number of iterations. The algorithm gets slow for a large number of
#'   objects. The speed can be improved by lowering the cooling parameter `"cool"`
#'   or increasing the minimum temperature `"tmin"`.
#'   However, this will decrease the seriation quality.
#'
#'    Directly minimizes the **linear seriation criterion (LS).**
#'
#' - **Complete Enumeration:** `"Enumerate"`
#'
#'    This method finds the optimal permutation given a seriation criterion by complete enumeration
#'    of all permutations.
#'    The criterion is specified as the `control` parameters `"criterion"`.
#'    Default is the weighted gradient measure. Use `"verbose = TRUE"` to see
#'    the progress.
#'
#'    Note: The number of permutations for \eqn{n} objects is \eqn{n!}.
#'    Complete enumeration is only possible for tiny problems (<10 objects) and is limited on most systems
#'    to a problem size of up to 12 objects.
#'
#' - **Gradient measure seriation by branch-and-bound:** `"BBURCG"`, `"BBWRCG"` (Brusco and Stahl 2005)
#'
#'    The method uses branch-and-bound to minimize the
#'    **unweighted gradient measure** (`"BBURCG"`) and the
#'    **weighted gradient measure** (`"BBWRCG"`).
#'    This type of optimization is only feasible for a small number of objects (< 50 objects).
#'
#'    For BBURCG, the control parameter `"eps"` can be used to relax the problem by defining
#'    that a distance needs to be eps larger to count as a violation. This relaxation will improve the speed,
#'    but miss some Robinson events. The default value is 0.
#'
#' - **Genetic Algorithm:** `"GA"`
#'
#'   The GA code has to be first registered. A detailed description can
#'   be found on the manual page for [register_GA()].
#'
#' - **Quadratic assignment problem seriation:**
#'    `"QAP_LS"`, `"QAP_2SUM"`, `"QAP_BAR"`, `"QAP_Inertia"` (Hahsler, 2017)
#'
#'   Formulates the seriation problem as a quadratic assignment problem and applies a
#'   simulated annealing solver to find a good solution.
#'   These methods minimize the
#'   **Linear Seriation Problem** (LS) formulation (Hubert and Schultz 1976),
#'   the **2-Sum Problem** formulation (Barnard, Pothen, and Simon 1993), the
#'   **banded anti-Robinson form** (BAR), or the **inertia criterion**.
#'
#'   `control` parameters are passed on to [qap::qap()].
#'   An important parameter is `rep` to return the best result from the
#'   given number of repetitions with random restarts. The default is 1, but bigger
#'   numbers result in better and more stable results.
#'
#' - **General Simulated Annealing:** `"GSA"`
#'
#'   Implement simulated annealing similar to the ARSA method. However, it
#'   can optimize
#'   for any criterion measure defined in \pkg{seriation}. By default, the
#'   algorithm optimizes for the raw gradient measure, and is warm started with the
#'   result of spectral seriation (2-Sum problem) since Hahsler (2017) shows that
#'   2-Sum solutions are similar to solutions for the gradient measure.
#'   Use `warmstart = "random"` for no warm start.
#'
#'   The initial temperature `t0` and minimum temperature `tmin` can be set. If
#'   `t0` is not set, then it is estimated by sampling uphill moves and setting
#'   `t0` such that the median uphill move have a probability
#'    of `tinitialaccept`.
#'    Using the cooling rate `cool`, the number of iterations
#'    to go for `t0` to `tmin` is calculated.
#'
#'   Several popular local neighborhood functions are
#'   provided, and new ones can be defined (see [LS]). Local moves are tried in each
#'   iteration `nlocal` times the number of objects.
#'
#'   Note that this is an R implementation repeatedly calling the criterion funciton
#'   which is very slow.
#'
#' - **Stochastic gradient descent:** `"SGD"`
#'
#'   Starts with a solution and then performs stochastic gradient descent to find
#'   a close-by local optimum given a specified criterion.
#'
#'   Important `control` parameters:
#'     - `"criterion"`: the criterion to optimize
#'     - `"init"`: initial seriation (an order or the name of a seriation method)
#'     - `"max_iter"`: number of trials
#'
#' - **Spectral seriation:** `"Spectral"`, `"Spectral_norm"`  (Ding and He, 2004)
#'
#'   Spectral seriation uses a relaxation to minimize the **2-Sum Problem**
#'   (Barnard, Pothen, and Simon, 1993). It uses the order of the Fiedler vector
#'   of the similarity matrix's (normalized) Laplacian.
#'
#'   Spectral seriation gives a good trade-off between seriation quality,
#'   and scalability (see Hahsler, 2017).
#'
#' - **Traveling salesperson problem solver:** `"TSP"`
#'
#'   Uses a traveling salesperson problem solver to minimize the
#'   **Hamiltonian path length**. The solvers in \pkg{TSP} are used (see
#'   [TSP::solve_TSP()]). The solver method can be passed on via the `control`
#'   argument, e.g., `control = list(method = "two_opt")`. Default is the est
#'   of 10 runs of arbitrary insertion heuristic with 2-opt improvement.
#'
#'   Since a tour returned by a TSP solver is a connected circle and we are
#'   looking for a path representing a linear order, we need to find the best
#'   cutting point.  Climer and Zhang (2006) suggest adding a dummy city with
#'   equal distance to each other city before generating the tour. The place of
#'   this dummy city in an optimal tour with minimal length is the best cutting
#'   point (it lies between the most distant cities).
#'
#'
#'  **Other Methods**
#'
#' - **Identity permutation:** `"Identity"
#'
#' - **Reverse Identity permutation:** `"Reverse"
#'
#' - **Random permutation:** `"Random"`
#'
#' - **Rank-two ellipse seriation:** `"R2E"`  (Chen 2002)
#'
#'   Rank-two ellipse seriation starts with generating a sequence of correlation matrices
#'   \eqn{R^1, R^2, \ldots}. \eqn{R^1} is the correlation matrix of the original
#'   distance matrix \eqn{D} (supplied to the function as `x`), and
#'   \deqn{R^{n+1} = \phi R^n,} where \eqn{\phi} calculates the correlation
#'   matrix.
#'
#'   The rank of the matrix \eqn{R^n} falls with increasing \eqn{n}. The first
#'   \eqn{R^n} in the sequence, which has a rank of 2 is found. Projecting all
#'   points in this matrix on the first two eigenvectors, all points fall on an
#'   ellipse. The order of the points on this ellipse is the resulting order.
#'
#'   The ellipse can be cut at the two interception points (top or bottom) of the
#'   vertical axis with the ellipse. In this implementation, the topmost cutting
#'   point is used.
#'
#' - **Sorting Points Into Neighborhoods:** `"SPIN_STS"`, `"SPIN_NH"` (Tsafrir, 2005)
#'
#'   Given a weight matrix \eqn{W}, the SPIN algorithms try to
#'   minimize the energy for a permutation (matrix \eqn{P}) given by \deqn{F(P) =
#'   tr(PDP^TW),} where \eqn{tr} denotes the matrix trace.
#'
#'   `"SPIN_STS"` implements the Side-to-Side algorithm, which tries to push
#'   out large distance values. The default weight matrix suggested in the paper
#'   with \eqn{W=XX^T} and \eqn{X_i=i-(n+1)/2} is used. We run the algorithm form
#'   `step` (25) iteration and restart the algorithm `nstart` (10) with
#'   random initial permutations (default values in parentheses).
#'
#'   `"SPIN_NH"` implements the neighborhood algorithm (concentrate low
#'   distance values around the diagonal) with a Gaussian weight matrix
#'   \eqn{W_{ij} = exp(-(i-j)^2/n\sigma)}, where \eqn{n} is the size of the
#'   dissimilarity matrix and \eqn{\sigma} is the variance around the diagonal
#'   that control the influence of global (large \eqn{\sigma}) or local (small
#'   \eqn{\sigma}) structure.
#'
#'   We use the heuristic suggested in the paper for the linear assignment
#'   problem. We do not terminate as indicated in the algorithm but run all the
#'   iterations since the heuristic does not guarantee that the energy is
#'   strictly decreasing. We also implement the heuristic "annealing" scheme
#'   where \eqn{\sigma} is successively reduced. The parameters in `control`
#'   are `sigma` which can be a single value or a decreasing sequence
#'   (default: 20 to 1 in 10 steps), and `step`, which defines how many update
#'   steps are performed before for each value of `alpha`. Via
#'   `W_function` a custom function to create \eqn{W} with the function
#'   signature `function(n, sigma, verbose)` can be specified.
#'
#' - **Visual Assessment of (Clustering) Tendency:** `"VAT"` (Bezdek and Hathaway, 2002).
#'
#'   Creates an order based on Prim's algorithm for finding a minimum spanning
#'   tree (MST) in a weighted connected graph representing the distance matrix.
#'   The order is given by the order in which the nodes (objects) are added to
#'   the MST.
#'
#'
#' **Seriation methods for matrices (matrix)**
#'
#' Two-mode two-way data are general matrices.
#' Some methods also require that the matrix is positive.
#' Data frames and contingency tables ([base::table])
#' are converted into a matrix. However, the
#' default methods are different.
#'
#' Some methods find the row and column order simultaneously,
#' while others calculate them independently.
#' Currently, the
#' following methods are implemented for `matrix`:
#'
#' **Seriating rows and columns simultaneously**
#'
#' Row and column order influence each other.
#'
#' - **Bond Energy Algorithm:** `"BEA"`  (McCormick, 1972).
#'
#'   The algorithm tries to maximize a non-negative matrix's
#'   **Measure of Effectiveness.**
#'   Due to the definition of this measure, the tasks of
#'   ordering rows and columns are separable and can be solved independently.
#'
#'   A row is arbitrarily placed; then, rows are positioned one by one. When this
#'   is completed, the columns are treated similarly. The overall procedure
#'   amounts to two approximate traveling salesperson problems (TSP), one on the
#'   rows and one on the columns. The so-called `best insertion strategy' is
#'   used: rows (or columns) are inserted into the current permuted list of rows
#'   (or columns). Several consecutive runs of the algorithm might improve the
#'   energy.
#'
#'   Arabie and Hubert (1990) question its use with non-binary data if
#'   the objective is to find a seriation or one-dimensional orderings of rows
#'   and columns.
#'
#'   Fionn Murtagh implemented the BEA code used in this package.
#'
#' - **TSP to optimize the Measure of Effectiveness**: `"BEA_TSP"` (Lenstra 1974).
#'
#'   Distances between rows are calculated for a \eqn{M \times N} data matrix as
#'   \eqn{d_{jk} = - \sum_{i=1}^{i=M} x_{ij}x_{ik}\ (j,k=0,1,...,N)}. Distances
#'   between columns are calculated the same way from the transposed data matrix.
#'
#'   Solving the two TSP using these distances optimizes the measure of
#'   effectiveness. BEA can be seen as a simple, suboptimal TSP method.
#'
#'   `control` parameter:
#'      - `"method"`: a TSP solver method (see [TSP::solve_TSP()]).
#'
#'
#' - **Correspondence analysis** `"CA"`
#'
#'   This function is designed to help simplify a mosaic plot or other displays of a
#'   matrix of frequencies.  It calculates a correspondence analysis of the matrix and
#'   an order for rows and columns according to the scores on a correspondence analysis dimension.
#'
#'   This is the default method for contingency tables.
#'
#'   `control` parameters:
#'     - `"dim"`: CA dimension used for reordering.
#'     - `"ca_param"`: List with parameters for the call to [ca::ca()].
#'
#'  **Seriating rows and columns separately using dissimilarities**
#'
#' - **Heatmap seriation:** `"Heatmap"`
#'
#'   Calculates distances between
#'   rows and between columns and then applies seriation so each. This is
#'   the default method for data frames.
#'
#'    `control` parameter:
#'      - `"seriation_method"`: a list with row and column seriation methods.
#'          The special method `"HC_Mean"` is available to use hierarchical clustering
#'          with reordering the leaves by the row/column means (see [stats::heatmap()]).
#'          Defaults to optimal leaf ordering `"OLO"`.
#'      -  `"seriation_control"`: a list with control parameters for row and column
#'          seriation methods.
#'      - `"dist_fun"`: specify the distance calculation as a function.
#'      - `"scale"`: `"none"`, `"row"`, or `"col"`.
#'
#'
#' **Seriate rows using the data matrix**
#'
#' These methods need access to the data matrix instead of dissimilarities to
#' reorder objects (rows). Columns can also be reorderd by applying the same technique
#' to the transposed data matrix.
#'
#' - **Order along the 1D locally linear embedding:** `"LLE"`
#'
#'  Performs 1D the non-linear dimensionality reduction method locally linear embedding
#'  (see [lle()]).
#'
#' - **Order along the first principal component:** `"PCA"`
#'
#'   Uses the projection of the data on its first principal component to
#'   determine the order of rows. Performs the same procedure on the transposed
#'   matrix to obtain the column order.
#'
#'   Note that for a distance matrix calculated from `x` with Euclidean
#'   distance, this method minimizes the least square criterion.
#'
#' - **Angular order of the first two PCA components:** `"PCA_angle"`
#'
#'   For rows, projects the data on the first two principal components
#'   and then orders by the angle in this space. The order is split by the larges
#'   gap between adjacent angles. A similar method was suggested by
#'   Friendly (2002) to order variables in correlation matrices
#'   by angles of first two eigenvectors. PCA also computes the eigenvectors
#'   of the covariance matrix of the data.
#'
#'   Performs the same process on the
#'   transposed matrix for the column order.
#'
#' **Other methods**
#'
#' - **Angular order of the first two eigenvectors:** `"AOE"` (Friendly 2002)
#'
#'    This method reordered correlation matrices by the angle in the space
#'    spanned by the two largest eigenvectors of the matrix. The order is split
#'    by the largest angle gap. This is the original method proposed by
#'    Friendly (2002).
#'
#' - **By row/column mean:** `"Mean"`
#'
#'    A transformation can be applied before calculating the means.
#'    The function is specified as control
#'    parameter `"transformation"`. Any function that takes as an input a
#'    matrix and returns the transformed matrix can be used. Examples
#'    are `scale` or `\(x) x^.5`.
#'
#'
#' - **Identity permutation:** `"Identity"`
#'
#' - **Reverse Identity permutation:** `"Reverse"`
#'
#' - **Random permutation:** `"Random"`
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
#' @param margin an integer vector giving the margin indices (dimensions) to be
#' seriated. For example, for a matrix, `1` indicates rows, `2`
#' indicates columns, `c(1 ,2)` means rows and columns.
#' Unseriated margins return the identity seriation order for that margin.
#' @param rep number of random restarts for randomized methods.
#'  Uses [seriate_rep()].
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
#' D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation
#' for application to visualization.
#' _Journal of Computational and Graphical Statistics,_ **24**(1), 1--25.
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
#' Kruskal, J.B. (1964). Nonmetric multidimensional scaling: a numerical method.
#' _Psychometrika_ **29**, 115--129.
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
#' Tenenbaum, J.B., de Silva, V. & Langford, J.C. (2000)
#' A global network framework for nonlinear dimensionality reduction.
#' _Science_ **290**, 2319-2323.
#'
#' Tsafrir, D., Tsafrir, I., Ein-Dor, L., Zuk, O., Notterman, D.A. and Domany,
#' E. (2005): Sorting points into neighborhoods (SPIN): data analysis and
#' visualization by ordering distance matrices, \emph{Bioinformatics,}
#' \bold{21}(10) 2301--8.
#' \doi{10.1093/bioinformatics/bti329}
#'
#' Sammon, J. W. (1969) A non-linear mapping for data structure analysis. _IEEE Trans. Comput._, **C-18** 401--409.
#' @keywords optimize cluster
#' @examples
#' # Show available seriation methods (for dist and matrix)
#' list_seriation_methods()
#'
#' # show the description for ARSA
#' get_seriation_method("dist", name = "ARSA")
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
#'
#' # plot in original order
#' pimage(corr, main = "Correlation matrix")
#'
#' # reorder the correlation matrix using the angle of eigenvectors
#' pimage(corr, order = "AOE", main = "Correlation matrix (AOE)")
#'
#' # we can also define a distance (we used d = sqrt(1 - r)) and
#' # then reorder the matrix (rows and columns) using any seriation method.
#' d <- as.dist(sqrt(1 - corr))
#' o <- seriate(d, method = "R2E")
#' corr_reordered <- permute(corr, order = c(o, o))
#' pimage(corr_reordered, main = "Correlation matrix (R2E)")
#' @export
seriate <- function(x, ...)
  UseMethod("seriate")

#' @export
seriate.default <- function(x, ...)
  stop(gettextf("seriate not implemented for class '%s'.",
    class(x)))
