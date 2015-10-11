\name{Zoo}
\alias{Zoo}
\docType{data}
\title{Zoo Data Set}
\description{
A database containing characteristics of different animals.
The database was created and donated by Richard S. Forsyth
and is available from the UCI Machine Learning Repository 
(Newman et al, 1998).
}
\usage{data("Zoo")}
\format{
  A data frame with 101 observations on the following 17 variables.
  \describe{
    \item{\code{hair}}{\{0, 1\}}
    \item{\code{feathers}}{\{0, 1\}}
    \item{\code{eggs}}{\{0, 1\}}
    \item{\code{milk}}{\{0, 1\}}
    \item{\code{airborne}}{\{0, 1\}}
    \item{\code{aquatic}}{\{0, 1\}}
    \item{\code{predator}}{\{0, 1\}}
    \item{\code{toothed}}{\{0, 1\}}
    \item{\code{backbone}}{\{0, 1\}}
    \item{\code{breathes}}{\{0, 1\}}
    \item{\code{venomous}}{\{0, 1\}}
    \item{\code{fins}}{\{0, 1\}}
    \item{\code{legs}}{Numeric (set of values: \{0, 2, 4, 5, 6, 8\})}
    \item{\code{tail}}{\{0, 1\}}
    \item{\code{domestic}}{\{0, 1\}}
    \item{\code{catsize}}{\{0, 1\}}
    \item{\code{class}}{a factor with levels \code{amphibian} \code{bird} \code{fish} \code{insect} \code{invertebrate} \code{mammal} \code{reptile}}
  }
}
%\details{
%}
\source{
D.J. Newman, S. Hettich, C.L. Blake and C.J. Merz (1998):
UCI Repository of machine learning databases,
\url{http://www.ics.uci.edu/~mlearn/MLRepository.html},
University of California, Irvine, Dept. of Information and Computer Sciences.
}
\examples{
data("Zoo")
x <- scale(Zoo[, -17])


d <- dist(x)
pimage(d)

order <- seriate(d, method = "tsp")
pimage(d, order)
}
\keyword{datasets}
