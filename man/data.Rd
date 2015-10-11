\name{seriation_data}
\alias{create_lines_data}
\alias{create_ordered_data}
\title{Create Simulated Data for Seriation Evaluation}
\description{
Several functions to create simulated data to evaluate different aspects
of seriation algorithms and criterion functions. 
}
\usage{
create_lines_data(n = 250)
create_ordered_data(n = 250, k = 2, size = NULL, spacing = 6,  path = "linear",
  sd1 = 1, sd2 = 0)
}
\arguments{
  \item{n}{number of data points to create.}
  \item{k}{number of Gaussian components.}
  \item{size}{relative size (number of points) of components (length of k). 
    If \code{NULL} then all components have the same size. }
  \item{spacing}{ space between the centers of components. The default of 6
    means that the components will barely touch at \code{ds1=1} 
    (3 standard deviations for each Gaussian component). }
  \item{path}{ Are the components arranged along a \code{"linear"} or 
  \code{"circular"} path? }
  \item{sd1}{ variation in the direction along the components. 
    A value greater than one means the
    components are mixing. }
  \item{sd2}{ variation perpendicular to the direction along the components. 
    A value greater than 0 will introduce anti-Robinson violation events. }
}
\details{
\code{create_lines_data} creates the lines data set used in 
for iVAT in Havens and Bezdeck (2012).


\code{create_ordered_data} is a versatile function which creates "orderable"
2D data using Gaussian components along a linear or circular path. The 
components are equally spaced (\code{spacing}) along the path. The 
default spacing of 6 ensures that 2 adjacent components with a standard 
deviation of one along the direction of the path will barely touch. The 
standard deviation along the path is set by \code{sd1}. The standard deviation
perpendicular to the path is set by \code{sd2}. A value larger than zero
will result in the data not being perfectly orderable (i.e., the
resulting distance matrix will not be a perfect pre-anti-Robinson matrix and
contain anti-Robinson violation events after seriation). Note that a circular
path always creates anti-Robinson violation since the circle has to be
broken at some point to create a linear order. 
}
\seealso{
\code{\link{seriate}},
\code{\link{criterion}},
\code{\link{VAT}}.
}
\references{
Havens, T.C. and Bezdek, J.C. (2012):
An Efficient Formulation of the Improved Visual Assessment of Cluster Tendency 
(iVAT) Algorithm, \emph{IEEE Transactions on Knowledge and Data Engineering,}
\bold{24}(5), 813--822. 
}
\author{Michael Hahsler}
\examples{
## lines data set from Havens and Bezdek (2011) 
x <- create_lines_data(250)
plot(x, xlim=c(-5,5), ylim=c(-3,3), cex=.2, col = attr(x, "id"))  
d <- dist(x)
pimage(d, seriate(d, "OLO_single"), col = bluered(100, bias=.5), key = TRUE)
  
## create_ordered_data can produce many types of "orderable" data

## perfect pre-Anti-Robinson matrix (with a single components)
x <- create_ordered_data(250, k = 1)  
plot(x, cex=.2, col = attr(x, "id"))  
d <- dist(x)
pimage(d, seriate(d, "MDS"), col = bluered(100, bias=.5), key = TRUE)  
  
## separated components
x <- create_ordered_data(250, k = 5)  
plot(x, cex=.2, col = attr(x, "id"))  
d <- dist(x)
pimage(d, seriate(d, "MDS"), col = bluered(100, bias=.5), key = TRUE)  
  
## overlapping components
x <- create_ordered_data(250, k = 5, sd1 = 2)  
plot(x, cex=.2, col = attr(x, "id"))
d <- dist(x)
pimage(d, seriate(d, "MDS"), col = bluered(100, bias=.5), key = TRUE)  
    
## introduce anti-Robinson violations (a non-zero y value)
x <- create_ordered_data(250, k = 5, sd1 = 2, sd2 = 5)  
plot(x, cex=.2, col = attr(x, "id"))  
d <- dist(x)
pimage(d, seriate(d, "MDS"), col = bluered(100, bias=.5), key = TRUE)  
  
## circular path (has always violations)  
x <- create_ordered_data(250, k = 5, path = "circular", sd1=2)  
plot(x, cex=.2, col = attr(x, "id"))
d <- dist(x)
pimage(d, seriate(d, "OLO"), col = bluered(100, bias=.5), key = TRUE)  

## circular path (with more violations violations)  
x <- create_ordered_data(250, k = 5, path = "circular", sd1=2, sd2=1)  
plot(x, cex=.2, col = attr(x, "id"))
d <- dist(x)
pimage(d, seriate(d, "OLO"), col = bluered(100, bias=.5), key = TRUE) 
}
\keyword{data}
