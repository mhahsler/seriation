\name{dissimilarity}
\alias{dist}
\alias{ser_dist}
\alias{ser_cor}
\alias{ser_align}
\title{Dissimilarities and Correlations Between Seriation Orders}
\description{
Calculates dissimilarities/correlations between seriation orders in a list.
}
\usage{
ser_cor(x, y = NULL, method = "spearman", reverse = FALSE, test = FALSE)
ser_dist(x, y = NULL, method = "spearman", reverse = FALSE)
ser_align(x, method = "spearman")
}
\arguments{
  \item{x}{set of seriation orders as a list with elements which can be 
  coerced into \code{ser_permutation_vector} objects. }
  \item{y}{ if not \code{NULL} then a single seriation order can be specified.   In this case \code{x} has to be a single seriation order and not a list. }
  \item{method}{ a character string with the name of the used measure. Available
    measures are:
    \code{"kendall"}, \code{"spearman"}, \code{"manhattan"}, 
    \code{"euclidean"}, \code{"hamming"}, and 
    \code{"ppc"} (positional proximity 
      coefficient).}
  \item{reverse}{a logical indicating if the orders should also be checked in 
    reverse order and the best value (highest correlation, lowed distance) is
    reported. This only affect ranking-based measures and not precedence
    invariant measures (e.g., ppc). } 
  \item{test}{a logical indicating if a correlation test should be performed. } 
}
\details{
\code{ser_cor} calculates the correlation between two sequences (orders).
Not that a seriation order and its reverse are identical and purely an artifact
due to the method that creates the order. This is a major difference to 
rankings.
For ranking-based correlation measures (Spearman and Kendall) 
the absolute value of the correlation is returned for \code{reverse = TRUE}
(in effect returning the correltation for the reversed order).
If \code{test = TRUE} then the appropriate test for association is performed
and a matrix with p-values is returned as the attribute \code{"p-value"}. Note
that no correction for multiple testing is performed.

For \code{ser_dist},
the correlation coefficients (Kendall's tau and Spearman's rho) are converted
into a dissimilarity by taking one minus the correlation value.
Note that Manhattan distance between the 
ranks in a linear order is equivalent to Spearman's footrule 
metric (Diaconis 1988). \code{reverse = TRUE} returns the pairwise minima
using also reversed orders.

The positional proximity coefficient (ppc) is a precedence invariant measure based on
the squared positional distances in two permutations (see Goulermas et al 2015).
We use the normalized value (i.e., the generalized correlation coefficient).
The similarity measure is converted into a dissimilarity via \eqn{1-ppc}.
For this precedence invariant measure \code{reverse} is ignored.

\code{ser_align} tries to normalize the direction in a list of seriations such 
that ranking-based methods can be used.
We add for each permutation also the reversed order to the set and then
use a modified version of Prim's
algorithm for finding a minimum spanning tree (MST) to choose if the original seriation order or its reverse should be used. We use the orders first added to 
the MST. Every time an order is added, its reverse is removed from the possible 
remaining orders.
}
\value{
\code{ser_dist} returns an object of class \code{dist}.
\code{ser_align} returns a new list with elements of class 
    \code{ser_permutation}.
}
\references{ 
P. Diaconis (1988): Group Representations in Probability and Statistics. Institute of Mathematical Statistics, Hayward, CA.

J.Y. Goulermas, A. Kostopoulos, and T. Mu (2015): A New Measure for Analyzing and Fusing Sequences of Objects. IEEE Transactions on Pattern Analysis and Machine Intelligence. Forthcomming.  
}
\seealso{
\code{\link{ser_permutation_vector}}
}
\author{Michael Hahsler}
\examples{
set.seed(1234)
## seriate dist of 50 flowers from the iris data set
data("iris")
x <- as.matrix(iris[-5])
x <- x[sample(1:nrow(x), 50),]
rownames(x) <- 1:50
d <- dist(x)

## Create a list of different seriations
methods <- c("HC_single", "HC_complete", "OLO", "GW", "R2E", "VAT", 
  "TSP", "Spectral", "SPIN", "MDS", "Identity", "Random")

os <- sapply(methods, function(m) {
  cat("Doing", m, "... ")
  tm <- system.time(o <- seriate(d, method = m))
  cat("took", tm[3],"s.\n")
  o
})

## Compare the methods using distances (default is based on 
## Spearman's rank correlation coefficient)
ds <- ser_dist(os)
hmap(ds, margin=c(7,7))

## Compare using actual correlation (reversed orders have low or 
## negative correlation!)
cs <- ser_cor(os)
hmap(cs, margin=c(7,7))

## Also check reversed seriation orders. 
## Now all but random and identity are highly positive correlated
cs2 <- ser_cor(os, reverse = TRUE)
hmap(cs2, margin=c(7,7))
  
## Use Manhattan distance of the ranks (i.e., Spearman's foot rule)
ds <- ser_dist(os, method="manhattan")
plot(hclust(ds))

## Also check reversed orders
ds <- ser_dist(os, method="manhattan", reverse = TRUE)
plot(hclust(ds))
}
\keyword{cluster}
