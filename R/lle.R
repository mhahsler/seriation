## lle is a simplified version from package lle
## by Holger Diedrich, Dr. Markus Abel

#' Locally Linear Embedding (LLE)
#'
#' Performs the non linear dimensionality reduction method locally linear embedding
#' proposed in Roweis and Saul (2000).
#'
#'
#' LLE tries to find a lower-dimensional projection which preserves distances
#' within local neighborhoods. This is done by (1) find for each object the
#' k nearest neighbors, (2) construct the LLE weight matrix
#' which represents each point as a linear combination of its neighborhood, and
#' (2) perform partial eigenvalue decomposition to find the embedding.
#'
#' The `reg` parameter allows the decision between different regularization methods.
#' As one step of the LLE algorithm, the inverse of the Gram-matrix \eqn{G\in R^{kxk}}
#' has to be calculated. The rank of \eqn{G} equals \eqn{m} which is mostly smaller
#' than \eqn{k} - this is why a regularization \eqn{G^{(i)}+r\cdot I} should be performed.
#' The calculation of regularization parameter \eqn{r} can be done using different methods:
#'
#' - `reg = 1`: standardized sum of eigenvalues of \eqn{G} (Roweis and Saul; 2000)
#' - `reg = 2` (default): trace of Gram-matrix divided by \eqn{k} (Grilli, 2007)
#' - `reg = 3`: constant value 3*10e-3
#'
#' @name lle
#' @aliases lle LLE
#'
#' @param x a matrix.
#' @param m dimensions of the desired embedding.
#' @param k number of neighbors.
#' @param reg regularization method. 1, 2 and 3, by default 2. See details.
#' @returns a matrix of vector with the embedding.
#' @author Michael Hahsler (based on code by Holger Diedrich and Markus Abel)
#' @references
#' Roweis, Sam T. and Saul, Lawrence K. (2000), Nonlinear Dimensionality
#' Reduction by Locally Linear Embedding,
#' _Science,_ **290**(5500), 2323--2326. \doi{10.1126/science.290.5500.2323}
#'
#' Grilli, Elisa (2007) Automated Local Linear Embedding with an application
#' to microarray data, Dissertation thesis, University of Bologna.
#' \doi{10.6092/unibo/amsdottorato/380}
#' @keywords cluster manip
#' @examples
#' data(iris)
#' x <- iris[, -5]
#'
#' # project iris on 2 dimensions
#' conf <- lle(x, m = 2, k = 30)
#' conf
#'
#' plot(conf, col = iris[, 5])
#'
#' # project iris onto a single dimension
#' conf <- lle(x, m = 1, k = 30)
#' conf
#'
#' plot_config(conf, col = iris[, 5], labels = FALSE)
#' @export
lle <-
  function(x,
           m,
           k,
           reg = 2) {
    nns <- find_nn_k(x, k)

    #calculate weights
    res_wgts <- find_weights(nns, x, m, reg)
    wgts <- res_wgts$wgts

    #compute coordinates
    y <- find_coords(wgts, nns, N = dim(x)[1], n = dim(x)[2], m)

    y
  }

find_coords <-
  function(wgts, nns, N, n, m)
  {
    W <- wgts
    M <- crossprod(diag(1, N) - W, diag(1, N) - W)
    eigen(M)$vectors[, c((N - m):(N - 1))] * sqrt(N)
  }


find_weights <-
  function(nns,
           x,
           m,
           reg = 2)
  {
    N <- dim(x)[1]
    n <- dim(x)[2]

    wgts <- 0 * matrix(0, N, N)
    #intrinsic dim
    intr_dim <- c()

    for (i in (1:N)) {
      #number of neighbours
      k <- sum(nns[i, ])

      #no neighbours (find_nn_k(k=0) or eps-neighbourhood)
      if (k == 0)
        next

      # calculate the differences  between xi and its neighbours
      Z <- matrix(c(t(x)) - c(t(x[i, ])), nrow = nrow(x), byrow = TRUE)
      Z <- matrix(Z[nns[i, ], ], ncol = n, nrow = k)

      #gram-matrix
      G <- Z %*% t(Z)

      #regularisation
      delta <- 0.1
      #calculate eigenvalues of G
      e <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
      #skip if all EV are null
      if (all(e == 0))
        next

      #choose regularisation method
      #see documentation
      if (reg == 1) {
        r <- delta * sum(utils::head(e, n - m)) / (n - m)
      } else if (reg == 2) {
        r <- delta ^ 2 / k * sum(diag(G))
      } else
        r <- 3 * 10 ^ -3

      #use regularization if more neighbors than dimensions!
      if (k > n)
        alpha <- r
      else
        alpha <- 0

      #regularization
      G <- G + alpha * diag(1, k)

      #calculate weights
      #using pseudoinverse ginv(A): works better for bad conditioned systems
      if (k >= 2)
        wgts[i, nns[i, ]] <- t(MASS::ginv(G) %*% rep(1, k))
      else
        wgts[i] <- G
      wgts[i, ] <- wgts[i, ] / sum(wgts[i, ])
    }


    return(list(
      x = x,
      wgts = wgts
    ))
  }

find_nn_k <-
  function(x, k) {
    nns <- as.matrix(dist(x))
    nns <- t(apply(nns, 1, rank))

    #choose the k+1 largest entries without the first (the data point itself)
    nns <= k + 1 & nns > 1
  }
