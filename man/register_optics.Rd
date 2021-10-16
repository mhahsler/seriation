\name{register_optics}
\alias{register_optics}
\alias{optics}
\alias{OPTICS}
\title{Register Seriation Based on OPTICS}
\description{
Use ordering points to identify the clustering structure (OPTICS) to create a seriation order.
}
\usage{
register_optics()
}
\details{
Registers the method \code{"optics"} for \code{seriate}. This method applies
the OPTICS ordering algorithm to create an ordering.

\bold{Note:} Package \pkg{dbscan} needs to be installed.
}
%\value{
%}
\references{
Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander (1999). OPTICS: Ordering Points To Identify the Clustering Structure. ACM SIGMOD international conference on Management of data. ACM Press. pp. 49-60. \doi{10.1145/304181.304187}
}
\seealso{
\code{\link{seriate}},
\code{\link[dbscan]{optics}} in \pkg{dbscan}.
}
\examples{
\dontrun{
register_optics()
get_seriation_method("dist", "optics")

d <- dist(random.robinson(50, pre=TRUE, noise=.1))

o <- seriate(d, method = "optics")
pimage(d, o)
}
}
\keyword{optimize}
\keyword{cluster}
