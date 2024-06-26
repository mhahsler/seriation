% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/register_DendSer.R
\name{register_DendSer}
\alias{register_DendSer}
\alias{DendSer}
\alias{dendser}
\title{Register Seriation Methods from Package DendSer}
\usage{
register_DendSer()
}
\value{
Nothing.
}
\description{
Register the DendSer dendrogram seriation method and the ARc criterion
(Earle and Hurley, 2015) for use with \code{\link[=seriate]{seriate()}}.
}
\details{
Registers the method \code{"DendSer"} for seriate. DendSer is a fast
heuristic for reordering dendrograms developed by Earle and Hurley (2015)
able to use different criteria.

\code{control} for \code{\link[=seriate]{seriate()}} with
method \code{"DendSer"} accepts the following parameters:
\itemize{
\item \code{"h"} or \code{"method"}: A dendrogram or a method for hierarchical clustering
(see \link{hclust}). Default: complete-link.
\item \code{"criterion"}: A seriation criterion to optimize (see
\code{list_criterion_methods("dist")}. Default: \code{"BAR"} (Banded
anti-Robinson from with 20\% band width).
\item \code{"verbose"}: a logical; print progress information?
\item \code{"DendSer_args"}: additional arguments for \code{\link[DendSer:DendSer]{DendSer::DendSer()}}.
}

For convenience, the following methods (for different cost functions) are
also provided:
\itemize{
\item \code{"DendSer_ARc"} (anti-robinson form),
\item \code{"DendSer_BAR"} (banded anti-Robinson form),
\item \code{"DendSer_LPL"} (lazy path length),
\item \code{"DendSer_PL"} (path length).
}

\strong{Note:} Package \pkg{DendSer} needs to be installed.
}
\examples{

\dontrun{
register_DendSer()
get_seriation_method("dist", "DendSer")

d <- dist(random.robinson(20, pre=TRUE))

## use Banded AR form with default clustering (complete-link)
o <- seriate(d, "DendSer_BAR")
pimage(d, o)

## use different hclust method (Ward) and AR as the cost function for
## dendrogram reordering
o <- seriate(d, "DendSer", control = list(method = "ward.D2", criterion = "AR"))
pimage(d, o)
}

}
\references{
D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation
for application to visualization. \emph{Journal of Computational and
Graphical Statistics,} \strong{24}(1), 1--25.
}
\seealso{
\code{\link[DendSer:DendSer]{DendSer::DendSer()}}

Other seriation: 
\code{\link{register_GA}()},
\code{\link{register_optics}()},
\code{\link{register_smacof}()},
\code{\link{register_tsne}()},
\code{\link{register_umap}()},
\code{\link{registry_for_seriaiton_methods}},
\code{\link{seriate}()},
\code{\link{seriate_best}()}
}
\author{
Michael Hahsler based on code by Catherine B. Hurley and Denise
Earle
}
\concept{seriation}
\keyword{cluster}
\keyword{optimize}
