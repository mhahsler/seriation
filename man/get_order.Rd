\name{get_order}
\alias{get_order}
\alias{get_order.integer}
\alias{get_order.hclust}
\alias{get_order.ser_permutation_vector}
\alias{get_order.ser_permutation}
\alias{get_rank}
\alias{get_permutation_matrix}
\title{Extracting Order Information from a Permutation Object}
\description{
Method to get the order information from an object of class
\code{ser_permutation} or \code{ser_permutation_vector}.
Order information can be extracted as an integer permutation vector,
a vector containing the object ranks or a permutation matrix.
}
\usage{
get_order(x, \ldots)
\method{get_order}{ser_permutation_vector}(x, \ldots)
\method{get_order}{ser_permutation}(x, dim = 1, \ldots)

get_rank(x, \ldots)

get_permutation_matrix(x, \ldots)
}
\arguments{
  \item{x}{ an object of class \code{ser_permutation} or
      \code{ser_permutation_vector}.}
  \item{dim}{ order information for which dimension should be returned?}
  \item{\ldots}{further arguments are ignored for \code{get_order}.
    For \code{get_rank} and
    for \code{get_permutation_matrix} the additional arguments are passed on to
    \code{get_order} (e.g., as \code{dim}).}
    }
\details{
\code{get_order} returns the seriation as an integer vector containing
  the order of the objects after permutation. That is, the index of the first,
  second, \eqn{..., n}-th object. These permuation vectors can directly be
  used to reorder objects using subsetting with \code{"["}.
  \emph{Note:} In \pkg{seriation} we usually use these
  order-based permutation vectors.

\code{get_rank} returns the seriation as an integer vector containing
  the rank/position for each objects in the permutation. That is,
  position of the first, second, \eqn{..., n}-th object after permutation.
  Note: Use \code{order()} to convert ranks back to an order.


\code{get_permutation_matrix} returns a \eqn{n \times n}{n x n} permutation matrix.
}
\value{
Returns an integer permutation vector/a permutation matrix.
}
\seealso{
\code{\link{ser_permutation_vector}},
\code{\link{ser_permutation}}
}
\author{Michael Hahsler}
\examples{
## permutation_vector
o <- ser_permutation_vector(sample(10))
o

get_order(o)
get_rank(o)
get_permutation_matrix(o)

## permutation
o2 <- ser_permutation(o, sample(5))
o2

get_order(o2, dim = 2)
get_rank(o2, dim = 2)
get_permutation_matrix(o2, dim = 2)
}
\keyword{manip}% at least one, from doc/KEYWORDS
