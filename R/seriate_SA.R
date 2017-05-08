#######################################################################
# seriation - Infrastructure for seriation
# Copyrigth (C) 2017 Michael Hahsler, Christian Buchta and Kurt Hornik
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


## Simulated annealing reimplimentation folowing arsa.f by Brusco et al.
## can use any criterion function


### neighborhood functions
LS_swap <- function(o, pos = sample.int(length(o), 2)) {
  tmp <- o[pos[1]]
  o[pos[1]] <- o[pos[2]]
  o[pos[2]] <- tmp
  o
}

### insert pos[1] in pos[2]
LS_insert <- function(o, pos = sample.int(length(o), 2)) {
  append(o[-pos[1]], o[pos[1]], after = pos[2]-1)
}

LS_reverse <- function(o, pos = sample.int(length(o), 2)) {
  o[pos[1]:pos[2]] <- o[pos[2]:pos[1]]
  o
}

LS_mixed <- function(o, pos = sample.int(length(o), 2)) {
  switch(sample.int(3, 1),
    LS_swap(o, pos),
    LS_insert(o, pos),
    LS_reverse(o, pos)
  )
}


seriate_sa <- function(x, control = NULL) {
  param <- .get_parameters(control, list(
    criterion = "Gradient_raw",
    init = "Spectral",    ## use "Random" for random init.
    localsearch = seriation::LS_insert,
    cool = 0.5,
    tmin = 0.0001,
    nlocal = 10,      ## try nlocal x n local search steps
    verbose = FALSE
  ))

  n <- attr(x, "Size")

  if(is.numeric(param$init)) {
    .check_dist_perm(x, order = param$init)
  }else{
    if(param$verbose) cat("\nObtaining initial solution via:",
      param$init, "\n")
    o <- get_order(seriate(x, method = param$init))
  }

  z <- criterion(x, o, method = param$criterion, force_loss = TRUE)
  if(param$verbose) cat("Initial z =", z,
    "(converted into loss if necessary)\n")

  iloop <- param$nlocal*n

  # find tmax (largest change for a move)
  znew <- replicate(iloop, expr = {
    criterion(x, param$localsearch(o), method = param$criterion,
      force_loss = TRUE)
  })

  tmax <- max(z-znew)
  if(tmax < 0) nloop <- 1L
  else nloop <- as.integer((log(param$tmin)-log(tmax))/log(param$cool))

  if(param$verbose) cat("Found tmax = ", tmax, "using", nloop, "iterations\n")

  zbest <- z
  temp <- tmax

  for(i in 1:nloop) {
    m <- 0L

    for(j in 1:iloop) {

      onew <- param$localsearch(o)
      znew <- criterion(x, onew, method = param$criterion, force_loss = TRUE)
      delta <- z-znew

      if(delta > 0 || runif(1) < exp(delta/temp)) {
        o <- onew
        z <- znew
        m <- m+1L
      }
    }

    if(param$verbose) {
      cat("temp = ", round(temp, 4), "\tz =", z,
        "\t performed moves = ", m, "/", iloop, "\n")
    }

    temp <- temp * param$cool
  }

  o
}

set_seriation_method("dist", "SA", seriate_sa,
  "Minimize a specified measure using simulated annealing (with warm start).")
