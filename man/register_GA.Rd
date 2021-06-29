\name{register_GA}
\alias{register_GA}
\alias{GA}
\alias{ga}
\alias{gaperm_mixedMutation}
\title{Register a Genetic Algorithm Seriation Method}
\description{
Register a GA-based seriation metaheuristic.
}
\usage{
register_GA()
}
\details{
Registers the method \code{"GA"} for \code{seriate}. This method
can be used to optimize any criterion in package \pkg{seriation}.
\code{control} for \code{seriate} with
method \code{"GA"} accepts the following parameters:
\describe{
\item{\code{"criterion"}}{ criterion to optimize. Default: BAR}
\item{\code{"suggestions"}}{ suggestions to warm start the GA. \code{NA}
  means no warm start. Default: TSP, QAP_LS and Spectral. }
\item{\code{"selection"}}{ Selection operator (see \pkg{GA}). Default: non-linear rank selection}
\item{\code{"crossover"}}{ Crossover operator (see \pkg{GA}). Default: ordered crossover (OX)}
\item{\code{"mutation"}}{ Mutation operator (see \pkg{GA}). Default:
  a mixture of the simple insertion (80\% chance) and simple
  inversion (20\% chance) operators.}
\item{\code{"pmutation"}}{ probaility for permutations. Default: .5}
\item{\code{"pcrossover"}}{ probability for crossover. Default: .2}
\item{\code{"popsize"}}{ the population size. Default: 100}
\item{\code{"maxiter"}}{ maximum number of generations. Default: 1000 }
\item{\code{"run"}}{ stop after \code{run} generations without improvement. Default: 50 }
\item{\code{"parallel"}}{ use multiple cores? Default: TRUE}
\item{\code{"verbose"}}{ Report progress? Default: TRUE}
}

The GA uses by default the ordered cross-over (OX) operator.
For mutation, the GA uses a mixture of simple insertion and simple inversion
operators. This mixed operator is created using
\code{seriation::gaperm_mixedMutation(ismProb = .8)}, where \code{ismProb} is
the probability that the simple insertion mutation operator is used.
See package \pkg{GA} for a description of other available
cross-over and mutation operators for permutations. The appropriate
operator functions in \pkg{GA} start with \code{gaperm_}.

We warm start the GA using \code{"suggestions"} given by several
heuristics. Set \code{"suggestions"} to \code{NA} to start with a purely
random initial population.

\bold{Note:} Package \pkg{GA} needs to be installed.
}
%\value{
%}
\references{
Luca Scrucca (2013): GA: A Package for Genetic Algorithms in R.
\emph{Journal of Statistical Software,} \bold{53}(4), 1--37.
URL \doi{10.18637/jss.v053.i04}.
}
\author{
Michael Hahsler
}
\seealso{
\code{\link{seriate}},
\code{\link[GA]{ga}} in \pkg{GA}.
}
\examples{
\dontrun{
register_GA()
list_seriation_methods("dist")

d <- dist(random.robinson(50, pre=TRUE, noise=.1))

## use default settings: Banded AR form
o <- seriate(d, "GA")
pimage(d, o)

## optimize for linear sertiation criterion (LS)
o <- seriate(d, "GA", control = list(criterion = "LS"))
pimage(d, o)

## no warm start
o <- seriate(d, "GA", control = list(criterion = "LS", suggestions = NA))
pimage(d, o)
}
}
\keyword{optimize}
\keyword{cluster}
