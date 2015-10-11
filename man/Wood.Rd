\name{Wood}
\alias{Wood}
\docType{data}
\title{Gene Expression Data for Wood Formation in Poplar Trees}
\description{
A data matrix containing a sample of the normalized gene expression data for
6 locations in the stem of Popla trees 
published in the study by Herzberg et al (2001).
The sample of 136 genes selected by 
Caraux and Pinloche (2005).
}
\usage{data(Wood)}
\format{
  The format is a 136 x 6 matrix.
}
\source{
The data was obtained from 
\url{http://www.atgc-montpellier.fr/permutmatrix/manual/Exemples/Wood/Wood.htm}.
}
\references{
Hertzberg M., H. Aspeborg, J. Schrader, A. Andersson, R.Erlandsson,
K. Blomqvist, R. Bhalerao, M. Uhlen, T. T. Teeri,
J. Lundeberg, Bjoern Sundberg, P. Nilsson and Goeran Sandberg (2001): 
A transcriptional roadmap to wood formation, 
\emph{PNAS,} \bold{98}(25), 14732--14737.

Caraux G. and Pinloche S. (2005):
PermutMatrix: a graphical environment to arrange gene expression profiles in optimal linear order, 
\emph{Bioinformatics,} \bold{21}(7) 1280--1281. 
}
\examples{
data(Wood)
head(Wood)
}
\keyword{datasets}
