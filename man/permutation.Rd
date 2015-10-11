\name{permutation}
\alias{ser_permutation}
\alias{permutation}
\alias{print.ser_permutation}
\alias{c.ser_permutation}
\alias{"[.ser_permutation"}
\title{Class ser_permutation -- A Collection of 
Permutation Vectors for Seriation}
\description{
The class \code{ser_permutation} is a collection of permutation vectors
(see class \code{ser_permutation_vector}), one for each dimension (mode) 
of the data to be permuted.
}
\usage{
## constructor
ser_permutation(x, \ldots)
}
\arguments{
  \item{x}{an object of class \code{ser_permutation_vector} or
      any object which can be converted into
      a object of class \code{ser_permutation} (e.g. an integer
      vector).}
  \item{\ldots}{permutation vectors for further dimensions}
    }
\details{
The basic functions \code{print}, \code{"["}, \code{"[["} and \code{c}
are provided.
}
\value{
An object of class \code{ser_permutation}.
}
\seealso{
\code{\link{ser_permutation_vector}},
\code{\link{get_order}},
\code{\link{get_permutation_matrix}}
}
\author{Michael Hahsler}
\examples{
o <- ser_permutation(1:5, 10:1)
o

## length (number of dimensions)
length(o)

## get permutation vector for 2nd dimension
get_order(o, 2)

## reverse dimensions
o[2:1]

## combine
o <- c(o, ser_permutation(1:15))
o

## get an individual permutation
o[[2]]

## reverse the order of a permutation
o[[2]] <- rev(o[[2]])
get_order(o,2)
}
\keyword{classes}% at least one, from doc/KEYWORDS
