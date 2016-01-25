#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2011 Michael Hahsler, Christian Buchta and Kurt Hornik
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## recognize Robinson structure


is.robinson <- function(x, anti = TRUE, pre = FALSE) {
  if(is.matrix(x) && !isSymmetric(unname(x)))
    stop("x needs to be a symmetric matrix!")

  d <- as.dist(x)
  if(!anti) d <- max(d) - d

  ## pre Robinson matrix can be perfectly seriated using
  ## spectral seriation!
  if(pre) d <- permute(d, seriate(d, method = "spectral"))

  unname(criterion(d, method = "AR_events") == 0)
}

random.robinson <- function(n, anti = TRUE, pre = FALSE, noise = 0) {

  if(noise < 0 | noise > 1) stop("noise has to be beween 0 and 1.")

  x <- runif(n)
  if(!pre) x <- sort(x)

  if(noise) x <- cbind(x, runif(n, min = 0, max = noise))

  m <- as.matrix(dist(x))

  if(!anti) m <- max(m)-m

  m
}

