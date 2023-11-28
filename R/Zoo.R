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
#' @source David Aha, Patrick Murphy, Christopher Merz, Eamonn Keogh,
#' Cathy Blake, Seth Hettich, David Newman, Arthur Asuncion, Moshe Lichman,
#' Dheeru Dua, Casey Graff (2023): UCI Machine Learning Repository,
#' \url{https://archive.ics.uci.edu/}, University of
#' California, Irvine.
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
