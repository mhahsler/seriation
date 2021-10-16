\name{register_umap}
\alias{register_umap}
\alias{umap}
\title{Register Seriation Based on 1D UMAP}
\description{
Use uniform manifold approximation and projection (UMAP) to embedd the data on the number line and create a seriation order.
}
\usage{
register_umap()
}
\details{
Registers the method \code{"umap"} for \code{seriate}. This method applies
1D UMAP to data represented by a distance matrix
and extracts the order from the 1D embedding.

\bold{Note:} Package \pkg{umap} needs to be installed.
}
%\value{
%}
\references{
McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
}
\seealso{
\code{\link{seriate}},
\code{\link[umap]{umap}} in \pkg{umap}.
}
\examples{
\dontrun{
register_umap()
get_seriation_method("dist", "umap")

d <- dist(random.robinson(50, pre=TRUE, noise=.1))

o <- seriate(d, method = "umap")
pimage(d, o)
}
}
\keyword{optimize}
\keyword{cluster}
