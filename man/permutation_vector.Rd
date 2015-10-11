\name{permutation_vector}
\alias{ser_permutation_vector}
\alias{permutation_vector}
\alias{print.ser_permutation_vector}
\alias{rev.ser_permutation_vector}
\alias{get_method}
\alias{length.ser_permutation_vector}
\title{Class ser_permutation_vector -- 
A Single Permutation Vector for Seriation}
\description{
The class \code{ser_permutation_vector} 
represents a single permutation vector. 
}
\usage{
## constructor
ser_permutation_vector(x, method = NULL)
}
\arguments{
  \item{x}{ an object which contains a permutation vector (currently an
      integer vector or an object of class \code{hclust}). The value \code{NA}
      creates an identity permutation. }
  \item{method}{ a string representing the method used to obtain the
      permutation vector}
    }
\details{
A permutation vector 
maps a set of \eqn{n} objects \eqn{\{O_1, O_2, ..., O_n\}}{{O_1, O_2, ..., O_n}} onto itself.
In \pkg{seriation} we represent a permutation \eqn{\pi}{\pi} 
as a vector which lists the objects in their permuted order. 
For example, the permutation vector \eqn{\langle3, 1, 2\rangle}{<3, 1, 2>} indicates that in first position is the object with index 3 then the object with index 1 and finally 
the object with index 2. 
A permutation vector can be extracted from a permuation vector object 
via \code{get_order()}. Such a permutation vector can be directly used
to subset the list of original objects with \code{"["} to apply the permutation.
\emph{Note:} An alternative way to specify a permutation is via a list of the ranks of the objects after permutation (see \code{get_rank()}).


\code{ser_permutation_vector} objects are usually packed into 
a \code{ser_permutation} object
which is a collection of \eqn{k} permutation vectors for \eqn{k}-mode data.

The constructor \code{ser_permutation_vector} 
checks if the permutation vector is valid 
(i.e. if all integers occur exactly once).

The following functions are implemented: \code{print}, \code{rev},
\code{length}, \code{get_order}, \code{get_rank}, \code{get_method}.
}
\value{
An object of class \code{ser_permutation_vector}.
}
\seealso{
\code{\link{ser_permutation}},
\code{\link{get_order}},
\code{\link{get_rank}},
\code{\link{get_permutation_matrix}},
\code{\link{permutation_vector2matrix}}.
}
\author{Michael Hahsler}
\examples{
p <- ser_permutation_vector(sample(10), "random")
p

## some methods
length(p)
get_method(p)
get_order(p)
get_rank(p)
get_permutation_matrix(p)

r <- rev(p)
r
get_order(r)

## create a indentity permutation vector (with unknown length)
ip <- ser_permutation_vector(NA)
ip 
}
\keyword{classes}% at least one, from doc/KEYWORDS
