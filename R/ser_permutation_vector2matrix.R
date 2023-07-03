#' Conversion Between Permutation Vector and Permutation Matrix
#'
#' Converts between permutation vectors and matrices.
#'
#' @family permutation
#'
#' @param x A permutation vector (any object that can be converted into a
#' permutation vector, e.g., a integer vector or a `hclust` object) or a
#' matrix representing a permutation. Arguments are checked.
#' @returns
#' - `permutation_vector2matrix()`: returns a permutation matrix.
#' - `permutation_matrix2vector()`: returns the permutation as a integer vector.
#'
#' @author Michael Hahsler
#' @keywords manip
#' @examples
#' ## create a random permutation vector
#' pv <- structure(sample(5), names = paste0("X", 1:5))
#' pv
#'
#' ## convert into a permutation matrix
#' pm <- permutation_vector2matrix(pv)
#' pm
#'
#' ## convert back
#' permutation_matrix2vector(pm)
#' @export
permutation_vector2matrix <- function(x) {
  x <- get_order(x)
  .valid_permutation_vector(x)

  n <- length(x)
  pm <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n)
    pm[i, x[i]] <- 1

  dimnames(pm) <- list(names(x), names(x))
  pm
}

#' @rdname permutation_vector2matrix
#' @export
permutation_matrix2vector <- function(x) {
  .valid_permutation_matrix(x)
  o <- apply(
    x,
    MARGIN = 1,
    FUN = function(r)
      which(r == 1)
  )
  o
}
