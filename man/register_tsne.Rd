% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/register_tsne.R
\name{register_tsne}
\alias{register_tsne}
\alias{tsne}
\alias{tSNE}
\title{Register Seriation Based on 1D t-SNE}
\usage{
register_tsne()
}
\value{
Nothing.
}
\description{
Use t-distributed stochastic neighbor embedding (t-SNE) for \code{\link[=seriate]{seriate()}}.
}
\details{
Registers the method \code{"tsne"} for \code{\link[=seriate]{seriate()}}. This method applies
1D t-SNE to a data matrix or a distance matrix and extracts the order
from the 1D embedding. To speed up the process, an initial embedding is
created using 1D multi-dimensional scaling (MDS) or principal
components analysis (PCA) which is improved by t-SNE.

The \code{control} parameter \code{"mds"} or \code{"pca"} controls if MDS (for distances)
or PCA (for data matrices) is used to create an
initial embedding. See \code{\link[Rtsne:Rtsne]{Rtsne::Rtsne()}} to learn about the other
available \code{control} parameters.

Perplexity is automatically set as the minimum between 30 and the number of
observations. It can be also specified using the control parameter
\code{"preplexity"}.

\strong{Note:} Package \pkg{Rtsne} needs to be installed.
}
\examples{

\dontrun{
register_tsne()

# distances
get_seriation_method("dist", "tsne")

data(SupremeCourt)
d <- as.dist(SupremeCourt)

o <- seriate(d, method = "tsne", verbose = TRUE)
pimage(d, o)

# look at the returned configuration and plot it
attr(o[[1]], "configuration")
plot_config(o)

# the t-SNE results are also available as an attribute (see ? Rtsne::Rtsne)
attr(o[[1]], "model")

## matrix
get_seriation_method("matrix", "tsne")

data("Zoo")
x <- Zoo

x[,"legs"] <- (x[,"legs"] > 0)

# t-SNE does not allow duplicates
x <- x[!duplicated(x), , drop = FALSE]

class <- x$class
label <- rownames(x)
x <- as.matrix(x[,-17])

o <- seriate(x, method = "tsne", eta = 10, verbose = TRUE)
pimage(x, o, prop = FALSE, row_labels = TRUE, col_labels = TRUE)

# look at the row embedding
plot_config(o[[1]], col = class)
}

}
\references{
van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing
High-Dimensional Data Using t-SNE. \emph{Journal of Machine Learning Research,}
\strong{9},
pp.2579-2605.
}
\seealso{
\code{\link[Rtsne:Rtsne]{Rtsne::Rtsne()}}

Other seriation: 
\code{\link{register_DendSer}()},
\code{\link{register_GA}()},
\code{\link{register_optics}()},
\code{\link{register_smacof}()},
\code{\link{register_umap}()},
\code{\link{registry_for_seriation_methods}},
\code{\link{seriate}()},
\code{\link{seriate_best}()}
}
\concept{seriation}
\keyword{cluster}
\keyword{optimize}
