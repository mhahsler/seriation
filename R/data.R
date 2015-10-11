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
    
### lines data set from iVAT paper

create_lines_data <- function(n=250) {
  n1 <- n/5*2
  n2 <- n/5
  n3 <- n/5*2
  
  x1 <- data.frame(x = runif(n1, -5, 5), y = rnorm(n1, mean = 2, sd = .1))
  x2 <- data.frame(x = runif(n2, -3, 3), y = rnorm(n2, mean = 0, sd = .1))
  x3 <- data.frame(x = runif(n3, -5, 5), y = rnorm(n3, mean = -2, sd = .1))
  id <- c(rep(1, times=n1), rep(2, times=n2), rep(3, times=n3))
  
  
  x <- rbind(x1, x2, x3)
  o <- sample(nrow(x))
  x <- x[o,]
  id <- id[o]
  
  rownames(x) <- 1:nrow(x)
  attr(x, "id") <- id
  
  x
}

### ordered data by Michael Hahsler (cite this package)
create_ordered_data <- function(n = 250, 
  k=2, size = NULL, spacing = 6, path = "linear",
  sd1=1, sd2=0) {
  
  if(k>n) stop("k needs to be less than n!")
  path <- match.arg(path, c("linear", "circular"))
  
  ## size
  if(is.null(size)) size <- rep(1, k)
  else if(length(size) != k) stop("length of size vector and k do not agree!")
  size <- round(size/sum(size) * n)
  size[1] <- n - sum(size[-1])
  
  ## create data
  ids <- rep(1:k, times = size)
  
  x <- data.frame(
    x = rnorm(n, mean = ids*spacing, sd = sd1),
    y = rnorm(n, mean = 0, sd = sd2)
  )

  ## transform
  if(path == "circular"){ 
    p <- k*spacing
    theta <- x[,1]/p * 2*pi
    r <-  p/(2*pi)+ x[,2]
    x <- cbind(x=r*sin(theta), y=r*cos(theta))
  }
  
  ## randomize order
  o <- sample(nrow(x))
  x <- x[o , , drop=FALSE]
  ids <- ids[o]
  attr(x, "id") <- ids
  
  x
}
