\name{VAT}
\alias{path_dist}
\alias{iVAT}
\alias{VAT}
\title{Visual Analysis for Cluster Tendency Assessment (VAT/iVAT)}
\description{
Implements Visual Analysis
for Cluster Tendency Assessment (VAT; Bezdek and Hathaway, 2002) and Improved
Visual Analysis for Cluster Tendency Assessment (iVAT; Wang et al, 2010).
}
\usage{
VAT(x, ...)
iVAT(x, ...)
path_dist(x)
}
\arguments{
  \item{x}{a \code{dist} object.}
  \item{...}{further arguments are passed on to \code{pimage}.}
}
\details{
\code{path_dist} redefines the distance between two objects as the minimum over
the largest distances in all possible paths between the objects as used for 
iVAT.
}
\seealso{
\code{\link{seriate}},
\code{\link{pimage}},
\code{\link{create_lines_data}}.
}
\references{
Bezdek, J.C. and Hathaway, R.J. (2002): VAT: a tool for visual assessment of (cluster) tendency. \emph{Proceedings of the 2002 International Joint 
Conference on Neural Networks (IJCNN '02)}, Volume: 3, 2225--2230. 
    
Havens, T.C. and Bezdek, J.C. (2012):
An Efficient Formulation of the Improved Visual Assessment of Cluster Tendency (iVAT) Algorithm, \emph{IEEE Transactions on Knowledge and Data Engineering,}   \bold{24}(5), 813--822. 

Wang L., U.T.V. Nguyen, J.C. Bezdek,
C.A. Leckie and K. Ramamohanarao (2010):
iVAT and aVAT: Enhanced Visual Analysis
for Cluster Tendency Assessment,
\emph{Proceedings of the PAKDD 2010, Part I, LNAI 6118,} 16--27.
}
\author{Michael Hahsler}
\examples{
## lines data set from Havens and Bezdek (2011) 
x <- create_lines_data(250)
plot(x, xlim=c(-5,5), ylim=c(-3,3), cex=.2)  
d <- dist(x)

## create regular VAT
VAT(d, colorkey = TRUE, main = "VAT")
## same as: pimage(d, seriate(d, "VAT"))  

## create iVAT which shows visually the three lines
iVAT(d, main = "iVAT")
## same as:
## d_path <- path_dist(d)
## pimage(d_path, seriate(d_path, "VAT"))  

## compare with dissplot (shows banded structures and relationship between 
## center line and the two outer lines)
dissplot(d, method="OLO_single", main = "Dissplot", col = bluered(100, bias = .5))

## compare with optimally reordered heatmap
hmap(d, method="OLO_single", main = "Heat map (opt. leaf ordering)", 
col = bluered(100, bias = .5))
}
\keyword{cluster}
\keyword{manip}
