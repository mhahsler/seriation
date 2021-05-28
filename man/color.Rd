\name{color_palettes}
\alias{greenred}
\alias{bluered}
\alias{grays}
\alias{greys}
\alias{reds}
\title{Different Useful Color Palettes}
\description{
Defines several color palettes for \code{pimage}, \code{dissplot}
and \code{hmap}.
}
\usage{
grays(n,  bias = 1, power = 1, ...)
greys(n,  bias = 1, power = 1, ...)
reds(n,  bias = 1, power = 1, ...)

bluered(n, bias = 1, power = 1, ...)
greenred(n,  bias = 1, power = 1, ...)
}
\arguments{
  \item{n}{ number of different colors produces. }
  \item{power}{ used to control how chroma and luminance is increased (1 = linear, 2 = quadratic, etc.)}
  \item{bias}{ a positive number. Higher values give more widely spaced colors at the high end. }
  \item{...}{ further parameters are passed on to
    \code{\link{sequential_hcl}} or \code{\link{diverging_hcl}}.}
}
\details{
The color palettes are created with
\code{sequential_hcl} and \code{diverging_hcl}
from package \pkg{colorspace}.

The two sequental palettes are: \code{reds} and
\code{grays} (or \code{greys}).

The two diverging palettes are:
\code{bluered} and \code{greenred}.
}
\value{
A vector with \code{n} colors.
}
\seealso{
  \code{\link[colorspace]{sequential_hcl}},
  \code{\link[colorspace]{diverging_hcl}},
  \code{\link{pimage}},
  \code{\link{dissplot}},
  \code{\link{hmap}}.
}
\author{Michael Hahsler}
\examples{
pimage(rbind(1:100))
pimage(rbind(1:100), col = greys(100, power = 2))
pimage(rbind(1:100), col = greys(100, bias = 2))
pimage(rbind(1:100), col = bluered(100))
pimage(rbind(1:100), col = bluered(100, power = .5))
pimage(rbind(1:100), col = bluered(100, bias = 2))
pimage(rbind(-50:100), col = greenred(20, bias = 1.5))

## choose your own color palettes
library(colorspace)
hcl_palettes(plot = TRUE)

## blues (with 20 shades)
pimage(rbind(1:100),
  col = colorspace::sequential_hcl(20, "Blues", rev = TRUE))
## blue to green
pimage(rbind(1:100),
  col = colorspace::diverging_hcl(100, "Tofino"))
}
\keyword{hplot}
