\name{SupremeCourt}
\alias{SupremeCourt}
\docType{data}
\title{Voting Patterns in the Second Rehnquist U.S. Supreme Court}
\description{
    Contains a (a subset of the) decisions for the stable 8-yr
    period 1995-2002 of the second Rehnquist Supreme Court.
    Decisions are aggregated to
    the joint probability for disagreement between judges.
}
\usage{data("SupremeCourt")}
\format{
  A square, symmetric 9-by-9 matrix with the joint probability for disagreement.
}
\references{
    Sirovich, L. (2003). A pattern analysis of the second Rehnquist
    U.S. Supreme Court. \emph{Proceedings of the National Academy of Sciences of the United States of America,} 100, 7432-7437. doi:10.1073/pnas.1132164100.
}
\author{Michael Hahsler}
\examples{
data("SupremeCourt")

SupremeCourt

d <- as.dist(SupremeCourt)
o <- seriate(d)

pimage(d, o)

plot(hclust(as.dist(SupremeCourt)))
}
\keyword{datasets}
