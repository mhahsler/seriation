\name{register_DendSer}
\alias{register_DendSer}
\alias{DendSer}
\alias{dendser}
\title{Register Seriation Methods from Package DendSer}
\description{
Register the DendSer dendrogram seriation method and
the ARc criterion (Earle and Hurley, 2015).
}
\usage{
register_DendSer()
}
\details{
Registers the method \code{"DendSer"} for \code{seriate}.
DendSer is a fast heuristic for reordering dendrograms developed by
Earle and Hurley (2015) able to use different criteria.
\code{control} for \code{seriate} with
method \code{"DendSer"} accepts the following parameters:
\describe{
\item{\code{"h"} or \code{"method"}}{A dendrogram or a method for
hierarchical clustering (see \code{hclust}). Default: complete-link.}
\item{\code{"cost"} or \code{"criterion"}}{Cost function to optimize.
Alternatively, you can specify a seriation
criterion (see \code{list_criterion_methods("dist")}).
Default: \code{"BAR"} (Banded anti-Robinson from with 20\% band width).}
\item{\code{"verbose"}}{print progress information.}
\item{\code{"DendSer_args"}}{additional arguments for \code{DendSer}.}
}

For convenience the following methods (for differnt cost functions)
are also provided:
\code{"DendSer_ARc"} (anti-robinson form),
\code{"DendSer_BAR"} (banded anti-Robinson form),
\code{"DendSer_LS"} (leaf seriation),
\code{"DendSer_PL"} (path length).

Note: Package \pkg{DendSer} needs to be installed.
}
%\value{
%}
\references{
D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation for application
to visualization. \emph{Journal of Computational and Graphical Statistics,}
\bold{24}(1), 1--25.
}
\author{
Michael Hahsler based on code by Catherine B. Hurley and Denise Earle
}
\seealso{
\code{\link{seriate}},
\code{\link[DendSer]{DendSer}} in \pkg{DendSer}.
}
\examples{
\dontrun{
register_DendSer()
list_seriation_methods("dist")

d <- dist(random.robinson(20, pre=TRUE))

## use Banded AR form with default clustering (complete-link)
o <- seriate(d, "DendSer_BAR")
pimage(d, o)

## use different hclust method (Ward) and AR as the cost function for
## dendrogram reordering
o <- seriate(d, "DendSer", control = list(method = "ward.D2", cost = "AR"))
pimage(d, o)
}
}
\keyword{optimize}
\keyword{cluster}
