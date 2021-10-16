\name{register_tsne}
\alias{register_tsne}
\alias{tsne}
\alias{tSNE}
\title{Register Seriation Based on 1D t-SNE}
\description{
Use t-distributed stochastic neighbor embedding (t-SNE) to create a seriation order.
}
\usage{
register_tsne()
}
\details{
Registers the method \code{"tsne"} for \code{seriate}. This method applies
1D t-SNE to data represented by a distance matrix
and extracts the order from the 1D embedding.
To speed up the process, an initial embedding is created using multi-dimensional scaling (MDS) which is
improved by t-SNE.

The \code{control} parameter \code{mds} controls if MDS is used to create an initial embedding.
See \code{\link[Rtsne]{Rtsne}} to learn about the other available
\code{control} parameters.

\bold{Note:} Package \pkg{Rtsne} needs to be installed.
}
%\value{
%}
\references{
van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing High-Dimensional Data Using t-SNE. Journal of Machine Learning Research, 9, pp.2579-2605.
}
\seealso{
\code{\link{seriate}},
\code{\link[Rtsne]{Rtsne}} in \pkg{Rtsne}.
}
\examples{
\dontrun{
register_tsne()
get_seriation_method("dist", "tsne")

d <- dist(random.robinson(50, pre=TRUE, noise=.1))

o <- seriate(d, method = "tsne")
pimage(d, o)
}
}
\keyword{optimize}
\keyword{cluster}
