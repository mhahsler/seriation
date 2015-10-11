\name{Robinson}
\alias{Robinson}
\alias{is.robinson}
\alias{random.robinson}
\title{Create and Recognize Robinson and Pre-Robinson Matrices}
\description{
Provides functions to create and 
recognize (anti) Robinson and pre-Robinson matrices.
A (anti) Robinson matrix has strictly decreasing (increasing) values when
moving away from the main diagonal. A pre-Robinson matrix is a matrix which
can be transformed into a perfect Robinson matrix using simultaneous 
permutations of rows and columns. 
}
\usage{
is.robinson(x, anti = TRUE, pre = FALSE)
random.robinson(n, anti = TRUE, pre = FALSE, noise = 0)
}
\arguments{
  \item{x}{a symmetric, positive matrix or a dissimilarity matrix 
  (a \code{dist} object). }
  \item{anti}{logical; check for anti Robinson structure? Note
  that for distances, anti Robinson structure is appropriate. }
  \item{pre}{logical; recognize/create pre-Robinson matrices. }
  \item{n}{ number of objects. } 
  \item{noise}{ zero means no noise. Noise more than zero results in 
    non-Robinson matrices. } 
}  
\details{
Note that the default matrices are anti Robinson matrices. This is done because
distance matrices (the default in R) are typically anti Robinson matrices
with values increasing when moving away from the diagonal.

Robinson matrices are recognized using the fact that they have zero anti Robinson events. For pre-Robinson matrices we use spectral seriation first since spectral seriation is guaranteed to perfectly reorder pre-Robinson matrices.

Random pre-Robinson matrices are generated 
using \eqn{dist(x)} where 
\eqn{x} is a vector of \code{n} values drawn from a random 
uniform distribution with range \eqn{[0,1]}. For Robinson matrices, 
\eqn{x} is sorted. 

For noise, randomly chosen \eqn{y} coordinates are added. The coordinates
are chosen from a normal distribution with a mean of zero and a standard deviation
given by the desired noise level. Reasonable levels are close to zero.
}
\value{
A single logical value.
}
%\seealso{
%}
\examples{
## create a perfect anti Robinson structure  
m <- random.robinson(10)
pimage(m)

is.robinson(m)

## permute the structure to make it not Robinsonian. However,
## it is still pre-Robinson.
o <- sample(10)
m2 <- permute(m, ser_permutation(o,o))
pimage(m2)

is.robinson(m2)
is.robinson(m2, pre = TRUE)

## create a binary random Robinson matrix (not anti Robinson)
m3 <- random.robinson(10, anti = FALSE) > .7 
pimage(m3)
is.robinson(m3, anti = FALSE)
}
%\keyword{manip}
