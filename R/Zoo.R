#' Zoo Data Set
#'
#' A database containing characteristics of different animals. The database was
#' created and donated by Richard S. Forsyth and is available from the UCI
#' Machine Learning Repository (Newman et al, 1998).
#'
#'
#' @name Zoo
#' @docType data
#' @format
#'   A data frame with 101 observations on the following 17 variables.
#' \describe{
#'   \item{\code{hair}}{a numeric vector}
#'   \item{\code{feathers}}{a numeric vector}
#'   \item{\code{eggs}}{a numeric vector}
#'   \item{\code{milk}}{a numeric vector}
#'   \item{\code{airborne}}{a numeric vector}
#'   \item{\code{aquatic}}{a numeric vector}
#'   \item{\code{predator}}{a numeric vector}
#'   \item{\code{toothed}}{a numeric vector}
#'   \item{\code{backbone}}{a numeric vector}
#'   \item{\code{breathes}}{a numeric vector}
#'   \item{\code{venomous}}{a numeric vector}
#'   \item{\code{fins}}{a numeric vector}
#'   \item{\code{legs}}{a numeric vector}
#'   \item{\code{tail}}{a numeric vector}
#'   \item{\code{domestic}}{a numeric vector}
#'   \item{\code{catsize}}{a numeric vector}
#'   \item{\code{class}}{a factor with levels \code{amphibian} \code{bird} \code{fish} \code{insect} \code{invertebrate} \code{mammal} \code{reptile}}
#' }
#' @source D.J. Newman, S. Hettich, C.L. Blake and C.J. Merz (1998): UCI
#' Repository of machine learning databases,
#' \url{https://www.ics.uci.edu/~mlearn/MLRepository.html}, University of
#' California, Irvine, Dept. of Information and Computer Sciences.
#' @keywords datasets
#' @examples
#'
#' data("Zoo")
#' x <- scale(Zoo[, -17])
#'
#'
#' d <- dist(x)
#' pimage(d)
#'
#' order <- seriate(d, method = "tsp")
#' pimage(d, order)
#'
NULL