\name{seriate}
\alias{seriate}
\alias{seriate.matrix}
\alias{seriate.dist}
\alias{seriate.array}
\alias{seriate.data.frame}
\alias{LS_insert}
\alias{LS_swap}
\alias{LS_reverse}
\alias{LS_mixed}
\title{Seriate Dissimilarity Matrices, Matrices or Arrays}
\description{
Tries to find an linear order for objects using data in form of a dissimilarity
matrix (two-way one mode data), a data matrix (two-way two-mode data) or a
data array (k-way k-mode data).
}
\usage{
\method{seriate}{dist}(x, method = "Spectral", control = NULL, \ldots)

\method{seriate}{array}(x, method = "PCA", control = NULL,
  margin = seq(length(dim(x))), \ldots)

\method{seriate}{matrix}(x, method = "PCA", control = NULL,
  margin = c(1,2), \ldots)

\method{seriate}{data.frame}(x, method = "Heatmap", control = NULL,
  margin = c(1,2), \ldots)
}
\arguments{
\item{x}{the data.}
\item{method}{ a character string with the name of the seriation method
(default: varies by data type).}
\item{control}{ a list of control options passed on to the seriation
algorithm.}
\item{margin}{ a vector giving the margin indices (dimensions) to be seriated. For example, for a matrix,
\code{1} indicates rows, \code{2} indicates columns, \code{c(1,2)}
indicates rows and columns.}
\item{\ldots}{ further arguments are added to the \code{control} list.}
}
\details{
Seriation methods are available via a registry.
See \code{\link{list_seriation_methods}} for help. Data frames are just a different reprentation of a
matrix and all seriation methods for matrix can be used for data frames. The default method for data frames is
heatmap seriation which calculates distances between rows and between columns and then applies seriation on these using
hierarchical clustering and optimal leaf ordering (OLO).

Many seriation methods (heuristically) optimize (minimize or maximize)
an objective function.
The value of the function for a given seriation can be calculated using
\code{\link{criterion}}. In this manual page, we only include the measure
which is optimized by each method using \bold{bold font}.
A definition of these measures can be found in the
\code{\link{criterion}} manual page.

\bold{Seriation methods by type of data:}

\bold{dist:} One-mode two-way data has to be provided as a dist object (not as a symmetric
matrix). Similarities have to be transformed into
dissimilarities. Currently, the following methods are implemented
(for a more detailed description and an experimental comparison see
\href{https://michael.hahsler.net/research/misc/EJOR_seriation_2016.pdf}{Hahsler (2017)}):

\describe{
\item{\code{"ARSA"}}{Anti-Robinson seriation by simulated annealing to minimize
the \bold{linear seriation criterion} (simulated annealing initialization used in Brusco et al 2008).
}

\item{\code{"BBURCG"}}{Anti-Robinson seriation by branch-and-bound to minimize
the \bold{unweighted gradient measure} (Brusco and Stahl 2005).
This is only feasible for a
relatively small number of objects.
}
\item{\code{"BBWRCG"}}{Anti-Robinson seriation by branch-and-bound to minimize
the \bold{weighted gradient measure} (Brusco and Stahl 2005). This is only feasible for a
relatively small number of objects.
}

\item{\code{"TSP"}}{Traveling salesperson problem solver to minimize
the \bold{Hamiltonian path length}.
The solvers in \pkg{TSP} are used (see \code{solve_TSP}).
The solver method can be passed on via the
\code{control} argument, e.g. \code{control = list(method = "two_opt")}.
Default is the est of 10 runs of
arbitrary insertion heuristic with 2-opt improvement.

Since a tour returned by a TSP solver is a connected circle and we are
looking for a path representing a linear order, we need to find the best
cutting point.  Climer and Zhang (2006) suggest to add a dummy city with
equal distance to each other city before generating the tour. The place
of this dummy city in an optimal tour with minimal length is the best
cutting point (it lies between the most distant cities).
}

\item{\code{"R2E"}}{Rank-two ellipse seriation (Chen 2002).

This method starts with generating a sequence of correlation matrices
\eqn{R^1, R^2, \ldots}. \eqn{R^1} is the correlation matrix
of the original distance matrix \eqn{D} (supplied to the function as
\code{x}),
and
\deqn{R^{n+1} = \phi R^n,}
where \eqn{\phi} calculates the
correlation matrix.

The rank of the matrix \eqn{R^n} falls with increasing \eqn{n}. The
first \eqn{R^n} in the sequence which has a rank of 2 is found.
Projecting all points in this matrix on the first two eigenvectors,
all points fall on an ellipse. The order of the points on this ellipse
is the resulting order.

The ellipse can be cut at the two interception points
(top or bottom) of the vertical axis with the ellipse.
In this implementation the top most cutting point is used.
}

\item{\code{"MDS"}, \code{"MDS_metric"}, \code{"MDS_nonmetric"},
\code{"MDS_angle"}}{Multidimensional scaling (MDS).

Use multidimensional scaling techniques to find an linear order
by minimizing \bold{stress}. Note
MDS algorithms used for a single dimension tend to end up in local optima
and unidimensional scaling (see Maier and De Leeuw, 2015)
would be more appropriate.
However, generally, ordering along the first component of
MDS provides good results.

By default, metric MDS (\code{cmdscale} in \pkg{stats}) is used.
In case of of general dissimilarities, non-metric MDS can be used.
The choices are \code{isoMDS} and \code{sammon} from \pkg{MASS}.
The method can be specified as the element \code{method}
(\code{"cmdscale"}, \code{"isoMDS"} or \code{"sammon"}) in \code{control}.

For convenience, \code{"MDS_metric"} performs \code{cmdscale} and
\code{"MDS_nonmetric"} performs \code{isoMDS}.

\code{"MDS_angle"} projects the data on the first two components
found by MDS and then orders by the angle in this space. The order
is split by the larges gap between adjacent angles. A similar method was
used for ordering correlation matrices by Friendly (2002).

}
\item{\code{"HC"}, \code{"HC_single"}, \code{"HC_complete"}, \code{"HC_average"},\code{"HC_ward"}}{Hierarchical clustering.

Using the order of the leaf nodes in a dendrogram obtained by hierarchical
clustering can be used as a very simple seriation technique.
This method applies hierarchical clustering (\code{hclust}) to \code{x}.
The clustering method can be given using a \code{"method"} element in
the \code{control} list. If omitted, the default \code{"average"} is
used.

For convenience the other methods are provided as shortcuts.
}


\item{\code{"GW"}, \code{"OLO"}}{Hierarchical
clustering (by default using average-link)
with additional leaf-node reordering to minimize
\bold{Hamiltonian path length (restricted)}.

A dendrogram (binary tree) has \eqn{2^{n-1}} internal nodes (subtrees) and
the same number of leaf orderings. That is, at each internal node the left
and right subtree (or leaves) can be swapped, or, in terms of a dendrogram,
be flipped.

Method \code{"GW"} uses an algorithm developed by Gruvaeus and Wainer
(1972) and implemented in package \pkg{gclus} (Hurley 2004).  The clusters are
ordered at each level so that the objects at the edge of each cluster are
adjacent to that object outside the cluster to which it is nearest. The method
produces an unique order.

Method \code{"OLO"} (Optimal leaf ordering, Bar-Joseph et al., 2001)
produces an optimal leaf ordering with respect to the
minimizing the sum of the distances along the (Hamiltonian) path connecting the
leaves in the given order. The time complexity of the algorithm is \eqn{O(n^3)}.
Note that non-finite distance values are not allowed.

Both methods start with a dendrogram created by \code{hclust}. As  the
\code{"method"} element in the \code{control} list a clustering method (default
\code{"average"}) can be specified. Alternatively, a \code{hclust} object can
be supplied using an element named \code{"hclust"}.

For convenience  \code{"GW_single"}, \code{"GW_average"},
\code{"GW_complete"}, \code{"GW_ward"} and
\code{"OLO_single"},    \code{"OLO_average"},
\code{"OLO_complete"},  \code{"OLO_ward"} are provided.
}

\item{\code{"VAT"}}{Visual Assessment of (Clustering) Tendency (Bezdek and Hathaway (2002)).

Creates an order based on Prim's algorithm for finding a minimum spanning
tree (MST) in a weighted connected graph representing the distance matrix.
The order is given by the order in which the nodes (objects) are added
to the MST.
}

\item{\code{"SA"}}{ Simulated Annealing for diverse criterion measures.

Implement simulated annealing similar to the ARSA method, however,
it works for any criterion measure defined in \pkg{seriation}.
By default the algorithm optimizes for raw gradient measure and is warm started with the result of spectral seriation (2-Sum problem)
since Hahsler (2017) shows that 2-Sum solutions are similar to
solutions for the gradient measure.

Local neighborhood functions are \code{LS_insert}, \code{LS_swap},
\code{LS_reverse}, and \code{LS_mix} (1/3 insertion, 1/3 swap and 1/3 reverse).
Any neighborhood function can be defined. It needs to take as the only argument
the order (integer vector) and return a random neighbor.

Note that this is an R implementation repeatedly calling criterion, and therefore is relatively slow.
}

\item{\code{"Spectral"}, \code{"Spectral_norm"}}{ Spectral seriation (Ding and He 2004).

Spectral seriation uses a relaxation to
minimize the \bold{2-Sum Problem} (Barnard, Pothen, and Simon 1993).
It uses the order of the Fiedler
vector of the similarity matrix's (normalized) Laplacian.

Spectral seriation gives a good trade-off between seriation quality, speed and
scalability (see Hahsler, 2017).
}

\item{\code{"SPIN_NH"}, \code{"SPIN_STS"}}{ Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005).
Given a weight matrix \eqn{W}, the algorithms try to minimize the energy for a permutation (matrix \eqn{P})
given by \deqn{F(P) = tr(PDP^TW),} where \eqn{tr} denotes the matrix
trace.

\code{"SPIN_STS"} implements the Side-to-Side algorithm which tries
to push out large distance values. The default weight matrix
suggested in the paper with \eqn{W=XX^T} and \eqn{X_i=i-(n+1)/2} is used.
We run the algorithm from \code{step} (25) iteration and restart the
algorithm \code{nstart} (10) with random initial permutations (default values in parentheses). Via \code{control} the parameters \code{step}, \code{nstart},
\code{X} and \code{verbose}.

\code{"SPIN_NH"} implements the neighborhood algorithm (concentrate
low distance values around the diagonal) with a
Gaussian weight matrix
\eqn{W_{ij} = exp(-(i-j)^2/n\sigma)}, where \eqn{n} is the size of the
dissimilarity matrix and \eqn{\sigma} is the variance around the diagonal
that control the influence
of global (large \eqn{\sigma}) or local (small \eqn{\sigma}) structure.

We use the heuristic suggested in the paper for the linear assignment problem.
We do not terminate as indicated in the algorithm, but run all the iterations since the heuristic does not guarantee that the energy is strictly decreasing.
We also implement the heuristic "annealing" scheme where \eqn{\sigma} is
successively reduced. The parameters in \code{control}
are \code{sigma} which can be a single value or a decreasing sequence
(default: 20 to 1 in 10 steps) and \code{step} which defines how many update
steps are performed before for each value of \code{alpha}.
Via \code{W_function} a custom function to create \eqn{W} with the function
signature \code{function(n, sigma, verbose)} can be specified.
The parameter \code{verbose} can be used to display progress information.
}

\item{\code{"QAP_LS"}, \code{"QAP_2SUM"}, \code{"QAP_BAR"}, \code{"QAP_Inertia"} }{ Quadratic assignment problem
formulations for seriation using a simulated annealing solver.
These methods minimize the \bold{Linear Seriation Problem} (LS)
formulation (Hubert and Schultz 1976),
the \bold{2-Sum Problem} formulation
(Barnard, Pothen, and Simon 1993),
the \bold{banded anti-Robinson form} (BAR) or
the \bold{inertia criterion.}

The parameters in \code{control} are passed on to \code{qap} in \pkg{qap}.
An important parameter is \code{rep} to return the best result out of
the given number of repetitions with random restarts. Default is 1, but bigger numbers result in better and more stable results.
}

\item{\code{"GA"}}{ Use a genetic algorithm to optimize for various criteria.
The GA code has to be first registered. A detailed description can be found in the manual page for \code{\link{register_GA}}.}

\item{\code{"DendSer"}}{ Use heuristic dendrogram seriation to optimize for various criteria.
The DendSer code has to be first registered. A detailed description can be found in the manual page for \code{\link{register_DendSer}}.}

\item{\code{"Identity"}}{ Produces an identity permutation. }

\item{\code{"Random"}}{ Produces a random permutation. }

}

\bold{matrix or data.frame:} Two-mode two-way data are general matrices.
Some methods also require that the matrix is positive.
Currently the following methods are implemented for matrix:
\describe{
\item{\code{"BEA"}}{Bond Energy Algorithm (BEA; McCormick 1972).
The algorithm tries to maximize the \bold{Measure of Effectiveness.}
of a non-negative matrix. Due to the definition
of this measure, the tasks of ordering rows and columns is separable
and can  be solved independently.

A row is arbitrarily placed; then rows are positioned one by one. When
this is completed, the columns are treated similarly. The overall
procedure amounts to two approximate traveling salesperson problems (TSP),
one on the rows and one on the columns. The so-called `best insertion
strategy' is used: rows (or columns) are inserted into the current
permuted list of rows (or columns). Several consecutive runs of the
algorithm might improve the energy.

Note that Arabie and Hubert (1990) question its use with non-binary data
if the objective is to find a seriation or one-dimensional orderings of
rows and columns.

The BEA code used in this package was implemented by Fionn Murtagh.

In \code{control} as element \code{"rep"} the number of runs can be
specified. The results of the best run will be returned.
}

\item{\code{"BEA_TSP"}}{Use a TSP to optimize the \bold{Measure of Effectiveness}
(Lenstra 1974).

In \code{control} as element \code{"method"} a TSP solver method can be
specified (see package \pkg{TSP}).
}

\item{\code{"Heatmap"}}{ Heatmap seriation calculates distances between rows and between columns and then applies seriation on these using hierarchical clustering and optimal leaf ordering (method \code{"OLO"} for \code{dist}). }

\item{\code{"PCA"}, \code{"PCA_angle"}}{ Principal component analysis.

Uses the projection of the data on its first principal component to
determine the order.

Note that for a distance matrix calculated from \code{x} with Euclidean
distance, this methods minimizes the least square criterion.

\code{"PCA_angle"} projects the data on the first two principal components
and then orders by the angle in this space. The order
is split by the larges gap between adjacent angles. A similar method was
used for ordering correlation matrices by Friendly (2002).
}

\item{\code{"Identity"}}{ Produces an identity permutation. }

\item{\code{"Random"}}{ Produces a random permutation. }

}

For array no built-in methods are currently available.

}
\value{
Returns an object of class \code{ser_permutation}.
}
\references{
Arabie, P. and L.J. Hubert (1990): The bond energy algorithm revisited,
\emph{IEEE Transactions on Systems, Man, and Cybernetics,}
\bold{20}(1), 268--274. \doi{10.1109/21.47829}

Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola. (2001):
Fast Optimal Leaf Ordering for Hierarchical Clustering.
\emph{Bioinformatics,} \bold{17}(1), 22--29. \doi{10.1093/bioinformatics/17.suppl_1.S22}

Barnard, S. T., A. Pothen, and H. D. Simon (1993): A Spectral Algorithm for Envelope Reduction of Sparse Matrices.
\emph{In Proceedings of the 1993
ACM/IEEE Conference on Supercomputing,} 493--502. Supercomputing '93.
New York, NY, USA: ACM. \url{https://ieeexplore.ieee.org/document/1263497}

Bezdek, J.C. and Hathaway, R.J. (2002): VAT: a tool for visual assessment of (cluster) tendency. \emph{Proceedings of the 2002 International Joint
Conference on Neural Networks (IJCNN '02)}, Volume: 3, 2225--2230.
\doi{10.1109/IJCNN.2002.1007487}

Brusco, M., Koehn, H.F., and Stahl, S. (2008): Heuristic Implementation of
Dynamic Programming for Matrix Permutation Problems in Combinatorial
Data Analysis. \emph{Psychometrika,} \bold{73}(3), 503--522. \doi{10.1007/s11336-007-9049-5}

Brusco, M., and Stahl, S. (2005):
\emph{Branch-and-Bound Applications in Combinatorial Data Analysis.}
New York: Springer. \doi{10.1007/0-387-28810-4}

Chen, C. H. (2002):  Generalized Association Plots: Information
Visualization via Iteratively Generated Correlation Matrices.
\emph{Statistica Sinica,} \bold{12}(1), 7--29.

Ding, C. and Xiaofeng He (2004): Linearized cluster assignment via spectral ordering. \emph{Proceedings of the Twenty-first International Conference on Machine learning (ICML '04)}. \doi{10.1145/1015330.1015407}

Climer, S. and Xiongnu Zhang (2006):
Rearrangement Clustering: Pitfalls, Remedies, and Applications,
\emph{Journal of Machine Learning Research,} \bold{7}(Jun),
919--943.

Friendly, M. (2002):
Corrgrams: Exploratory Displays for Correlation Matrices.
\emph{The American Statistician}, \bold{56}(4), 316--324.
\doi{10.1198/000313002533}

Gruvaeus, G. and Wainer, H. (1972):
Two Additions to Hierarchical Cluster Analysis,
\emph{British Journal of Mathematical and Statistical Psychology,}
\bold{25}, 200--206. \doi{10.1111/j.2044-8317.1972.tb00491.x}

Hahsler, M. (2017): An experimental comparison of seriation methods for one-mode two-way data. \emph{European Journal of Operational Research,} \bold{257}, 133--143. \doi{10.1016/j.ejor.2016.08.066}

Hubert, Lawrence, and James Schultz (1976):
Quadratic Assignment as a General Data Analysis Strategy.
\emph{British Journal of Mathematical and Statistical Psychology}
\bold{29}(2). Blackwell Publishing Ltd. 190--241.
\doi{10.1111/j.2044-8317.1976.tb00714.x}

Hurley, Catherine B. (2004):
Clustering Visualizations of Multidimensional Data.
\emph{Journal of Computational and Graphical Statistics,}
\bold{13}(4), 788--806.
\doi{10.1198/106186004X12425}

Lenstra, J.K (1974): Clustering a Data Array and the Traveling-Salesman
Problem, \emph{Operations Research,} \bold{22}(2) 413--414.
\doi{10.1287/opre.22.2.413}

Mair P., De Leeuw J. (2015). Unidimensional scaling.
In \emph{Wiley StatsRef: Statistics Reference Online,} Wiley, New York.
\doi{10.1002/9781118445112.stat06462.pub2}

McCormick, W.T., P.J. Schweitzer and T.W. White (1972):
Problem decomposition and data reorganization by a clustering technique,
\emph{Operations Research,} \bold{20}(5), 993--1009.
\doi{10.1287/opre.20.5.993}

Tsafrir, D., Tsafrir, I., Ein-Dor, L., Zuk, O., Notterman, D.A. and Domany, E. (2005): Sorting points into neighborhoods (SPIN): data analysis and visualization by ordering distance matrices, \emph{Bioinformatics,} \bold{21}(10) 2301--8.
\doi{10.1093/bioinformatics/bti329}
}
\seealso{
\code{\link{list_seriation_methods}},
\code{\link{criterion}},
\code{\link{register_GA}},
\code{\link{register_DendSer}},
\code{\link[TSP]{solve_TSP}} in \pkg{TSP},
\code{\link{hclust}} in \pkg{stats}.
}
\author{Michael Hahsler}
\examples{
# Show available seriation methods (for dist and matrix)
list_seriation_methods()

# Seriate as distance matrix (for the iris dataset)
data("iris")
x <- as.matrix(iris[-5])
x <- x[sample(1:nrow(x)), ]
d <- dist(x)

# Default seriation
order <- seriate(d)
order

pimage(d, main = "Random")
pimage(d, order, main = "Reordered")

# Compare seriation quality
rbind(
        random = criterion(d),
        reordered = criterion(d, order)
     )

# Seriate matrix
data("iris")
x <- as.matrix(iris[-5])

# To make the variables comparable, we scale the data
x <- scale(x, center = FALSE)

# The iris flowers are ordered by species in the data set
pimage(x, main = "original data", prop = FALSE)
criterion(x)

# Apply some methods
order <- seriate(x, method = "BEA_TSP")
pimage(x, order, main = "TSP to optimize ME", prop = FALSE)
criterion(x, order)

order <- seriate(x, method = "PCA")
pimage(x, order, main = "First principal component", prop = FALSE)
criterion(x, order)

order <- seriate(x, method = "heatmap")
pimage(x, order, main = "Heatmap seriation", prop = FALSE)
criterion(x, order)

# create a heatmap seriation manually by calculating
# distances between rows and between columns
order <- c(
    seriate(dist(x), method = "OLO"),
    seriate(dist(t(x)), method = "OLO")
)
pimage(x, order, main = "Heatmap seriation", prop = FALSE)
criterion(x, order)
}
\keyword{optimize}
\keyword{cluster}

