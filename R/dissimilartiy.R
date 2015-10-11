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


.dist_methods <- c("spearman", "kendall", "manhattan", "euclidean", "hamming",
  "ppc")

ser_cor <- function(x, y = NULL, method = "spearman", 
  reverse = FALSE, test=FALSE) { 
  ## Note: not all .dist_methods are implemented!
  method <- match.arg(tolower(method), .dist_methods)
  
  ## make sure everything is a permutation vector
  if(!is.null(y)) 
    x <- list(ser_permutation_vector(x), ser_permutation_vector(y))
  else x <- lapply(x, ser_permutation_vector)
  
  m <- .lget_rank(x)
  
  if(method == "ppc") {
    if(test) stop("No test for association available for PPC!")
    return(.ppc(x))
  }
    
  ## cor based methods 
  co <- cor(m, method = method)
  if(reverse) co <- abs(co) 

  ## add a test?
  if(test) {
    p <- outer(1:ncol(m), 1:ncol(m), FUN = 
      Vectorize(
        function(i, j) cor.test(m[,i], m[,j], method = method)$p.value))
    dimnames(p) <- dimnames(co)
    attr(co, "p-value") <- p
  }
  
  co
}

ser_dist <- function(x, y = NULL, method = "spearman", reverse = FALSE) {
  
  method <- match.arg(tolower(method), .dist_methods)
  
  ## make sure everything is a permutation vector
  if(!is.null(y)) 
    x <- list(ser_permutation_vector(x), ser_permutation_vector(y))
  else x <- lapply(x, ser_permutation_vector)
  
  if(!reverse) switch(method,
    spearman = as.dist(1-ser_cor(x, method="spearman", reverse = FALSE)),
    kendall = as.dist(1-ser_cor(x, method="kendal", reverse = FALSE)),
    
    ### Manhattan == Spearman's footrule  
    manhattan = dist(t(.lget_rank(x)), method="manhattan"),
    euclidean = dist(t(.lget_rank(x)), method="euclidean"),
    hamming   = .dist_hamming(t(.lget_rank(x))), 
    ppc = as.dist(1-ser_cor(x, method="ppc", reverse = FALSE))
  )
  
  else switch(method,
    spearman = as.dist(1-ser_cor(x, method="spearman", reverse = TRUE)),
    kendall =  as.dist(1-ser_cor(x, method="kendal", reverse = TRUE)),
    
    ### Manhattan == Spearman's footrule  
    manhattan = .find_best(dist(t(.lget_rank(.add_rev(x))), 
      method="manhattan")),
    euclidean = .find_best(dist(t(.lget_rank(.add_rev(x))), 
      method="euclidean")),
    hamming   = .find_best(.dist_hamming(t(.lget_rank(.add_rev(x))))),
    
    ### positional proximity coefficient is direction invariant
    ppc = as.dist(1-ser_cor(x, method="ppc", reverse = FALSE))
  )
}

ser_align <- function(x, method = "spearman") {
    if(!is.list(x) || any(!sapply(x, is, "ser_permutation_vector"))) 
      stop("x needs to be a list with elements of type 'ser_permutation_vector'")  
  
    .do_rev(x, .alignment(x, method=method))
}

.dist_hamming <- function(x) {
  n <- nrow(x)
  m <- matrix(nrow=n, ncol=n)
  for(i in seq_len(n))
    for(j in seq(i, n))
      m[j, i] <- m[i, j] <- sum(x[i,] != x[j,])
  mode(m) <- "numeric"
  dimnames(m) <- list(rownames(x), rownames(x))
  as.dist(m)
}

### make a permutation list into a rank matrix (cols are permutations)
.lget_rank <- function(x) sapply(x, get_rank)

### add reversed permutations to a list of permutations
.add_rev <- function(x) {
  os <- append(x, lapply(x, rev))
  names(os) <- c(labels(x), paste(labels(x), "_rev", sep=""))
  os
}

### reverses permutations in the list given a logical indicator vector
.do_rev <- function(x, rev) {
  for(i in which(rev)) x[[i]] <- rev(x[[i]])
  x
}

### finds the smallest distance in lists with reversed orders present 
.find_best <- function(d) {
  ### find smallest values
  m <- as.matrix(d)
  n <- nrow(m)/2
  m1 <- m[1:n, 1:n]
  m2 <- m[(n+1):(2*n), (n+1):(2*n)]
  m3 <- m[1:n, (n+1):(2*n)]
  m4 <- m[(n+1):(2*n), 1:n]
  as.dist(pmin(m1, m2, m3, m4))
}

### find largest values in matrix
.find_best_max <- function(d) {
  m <- as.matrix(d)
  n <- nrow(m)/2
  m1 <- m[1:n, 1:n]
  m2 <- m[(n+1):(2*n), (n+1):(2*n)]
  m3 <- m[1:n, (n+1):(2*n)]
  m4 <- m[(n+1):(2*n), 1:n]
  pmax(m1, m2, m3, m4)
}


### returns TRUE for sequences which should be reversed
.alignment <- function(x, method = "spearman") {
    if(!is.list(x) || any(!sapply(x, is, "ser_permutation_vector"))) stop("x needs to be a list with elements of type 'ser_permutation_vector'")  
  
    method <- match.arg(tolower(method), .dist_methods) 
    
    n <- length(x)
  
    ## calculate dist (orders + reversed orders)  
    d <- as.matrix(ser_dist(.add_rev(x), method=method, reverse=FALSE))
    diag(d) <- NA
    for(i in 1:n) { 
      d[i, n+i] <- NA
      d[n+i, i] <- NA  
    }
     
    ## start with closest pair
    take <- which(d == min(d, na.rm = TRUE), arr.ind = TRUE)[1,]
    #d[, c(take, (take+n) %% (2*n))] <- NA
    
    ## mark order and complement as taken
    d[, c(take, (take+n) %% (2*n))] <- Inf

    ## keep adding the closest 
    while(length(take) < n) {
      t2 <- which(d[take,] == min(d[take,], na.rm = TRUE), arr.ind = TRUE)[1, 2]
      #d[, c(t2,  (t2+n) %% (2*n))] <- NA
      
      ### closest to all
      #t2 <- which.min(colSums(d[take,], na.rm = T))
      d[, c(t2,  (t2+n) %% (2*n))] <- Inf
      
      take <- append(take, t2)
    }
  
    ## create indicator vector for the orders which need to be reversed
    take_ind <- logical(n)
    take_ind[take[take>n]-n] <- TRUE
    names(take_ind) <- names(x)
    take_ind
}
   
## Propositional Proximity Coefficient (1 - generalized corr. coef.)
## Goulermas, Kostopoulos and Mu, A new measure for analyzing and fusing 
## sequences of objects, IEEE Transactions on Pattern Analysis and Machine
## Intelligence, forthcomming.
##
## x,y ... permutation vectors (ranks)
.ppc_int <- function(x, y) {
  x <- get_rank(x)
  y <- get_rank(y)
  n <- length(x)
  
  #sum <- 0
  #for(j in 2:n) for(i in 1:(j-1)) sum <- sum + (x[i]-x[j])^2 * (y[i]-y[j])^2

  ## use fast matrix algebra instead
  Ax <- (x %*% rbind(rep_len(1, n)) - tcrossprod(cbind(rep_len(1, n)), x))^2
  Ay <- (y %*% rbind(rep_len(1, n)) - tcrossprod(cbind(rep_len(1, n)), y))^2
  ## note: Ay is symetric
  sum <- sum(diag(Ax %*% Ay)) 
  
  ## scale by theoretical maximum
  sum / (n^6/15 - n^4/6 + n^2/10) 
}

.vppc <- Vectorize(.ppc_int)
.ppc <- function(x) outer(x, x, .vppc)