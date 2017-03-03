\name{uniscale}
\alias{uniscale}
\alias{orderplot}
\title{Unidimensional Scaling from Seriation Results}
\description{
Performs (approximate) unidimensional scaling by first performin seriation
to obtain a permutation and the using the permutation to calulate the
configuration.
}
\usage{
uniscale(d, order = NULL, method = "QAP_LS", rep = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{a dissimilarity matrix.}
  \item{order}{a precomputed permutation order obtained via seriation.
  If \code{NULL}, then seriation is performed using the method
  specified in \code{method}.
}
  \item{method}{seriation method used if \code{o} is \code{NULL}.
}
  \item{rep}{Number of repetitions of the seriation heuristic.}
  \item{\dots}{additional arguments are passed on to the seriation method.}
}
\details{
Uses the method describes in Maier and De Leeuw (2015) to calculate the
minimum stress configuration for a, either a given or a
permutation created via a seriation method.

The code is similar to \code{uniscale} in \pkg{smacof}, but scales to larger
datasets since it does not check all permutations.
}
\value{
A vector with the fitted configuration.
}
\references{
Mair P., De Leeuw J. (2015). Unidimensional scaling. In Wiley StatsRef: Statistics Reference Online, Wiley, New York.
}
\author{
Michael Hahsler with code from Patrick Mair (from \pkg{smacof}).
}
\seealso{
\code{\link[smacof]{uniscale}} in \pkg{smacof}.
}
\examples{
data(Psych24)

# convert the correlation matrix into a distance matrix
d <- as.dist(1 - abs(Psych24))

sc <- uniscale(d)
sc

orderplot(sc)
}
\keyword{optimize}% use one of  RShowDoc("KEYWORDS")
