#' 2D Data Sets used for the CHAMELEON Clustering Algorithm
#'
#' Several 2D data sets created to evaluate the CHAMELEON clustering algorithm in
#' the paper by Karypis et al (1999).
#'
#' @name Chameleon
#' @aliases Chameleon chameleon chameleon_ds4 chameleon_ds5 chameleon_ds7
#' chameleon_ds8
#' @docType data
#' @family data
#' @format
#' `chameleon_ds4`: The format is a 8,000 x 2 data.frame.
#'
#' `chameleon_ds5`: The format is a 8,000 x 2 data.frame.
#'
#' `chameleon_ds7`: The format is a 10,000 x 2 data.frame.
#'
#' `chameleon_ds8`: The format is a 8,000 x 2 data.frame.
#' @references Karypis, G., EH. Han, V. Kumar (1999): CHAMELEON: A Hierarchical
#' Clustering Algorithm Using Dynamic Modeling, _IEEE Computer,_
#' **32**(8): 68--75.
#' \doi{10.1109/2.781637}
#' @keywords datasets
#' @examples
#' data(Chameleon)
#'
#' plot(chameleon_ds4, cex = .1)
#' plot(chameleon_ds5, cex = .1)
#' plot(chameleon_ds7, cex = .1)
#' plot(chameleon_ds8, cex = .1)
NULL

# link does not work
# @source The data was obtained from
# \url{http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download}
