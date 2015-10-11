\name{color_palettes}
\alias{greenred}
\alias{bluered}
\alias{grays}
\alias{greys}
\title{Different Useful Color Palettes}
\description{
Defines several color palettes for \code{pimage}, \code{dissplot}
and \code{hmap}.
}
\usage{
bluered(n, bias = 1)
greenred(n, bias = 1)
grays(n, power = 1)
greys(n, power = 1)
}
\arguments{
  \item{n}{ number of different colors produces. }
  \item{bias}{a positive number. Higher values give more widely 
    spaced colors at the high end.}
  \item{power}{ control parameter determining how luminance should be 
    increased (1 = linear, 2 = quadratic, etc.).}
}
\details{
\code{bluered} creates a blue-red color palette.

\code{greenred} creates a green-red color palette.

\code{greys} and \code{grays} creates gray scales.

See \code{colorRampPalette} to create your own color palettes.
}
\value{
A vector with \code{n} colors.
}
\seealso{
  \code{\link[grDevices]{colorRampPalette}},
  \code{\link{pimage}},
  \code{\link{dissplot}},
  \code{\link{hmap}}.
}
\author{Michael Hahsler}
\examples{
pimage(rbind(1:100))
pimage(rbind(1:100), col = greys(100, power=2))
pimage(rbind(1:100), col = bluered(100))
pimage(rbind(1:100), col = bluered(100, bias = 2))
pimage(rbind(-100:100), col = greenred(10))

## create your own color palettes 
## red to green (with 10 colors)
pimage(rbind(1:100), 
  col = colorRampPalette(colors = c("red", "yellow", "green"))(10))
## white to blue (with 100 colors)
pimage(rbind(1:100), 
  col = colorRampPalette(colors = c("white", "blue"))(100))
}
\keyword{hplot}
