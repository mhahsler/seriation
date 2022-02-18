#' 2D Data Sets used for the CHAMELEON Clustering Algorithm
#'
#' Several 2D data sets used to evaluate the CHAMELEON clustering algorithm in
#' the paper by Karypis et al (1999) and used by iVAT, an ordering-based tool
#' to asses cluster tendency (Havens and Bezdek, 2012).
#'
#' @name Chameleon
#' @aliases Chameleon chameleon chameleon_ds4 chameleon_ds5 chameleon_ds7
#' chameleon_ds8
#' @docType data
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
#'
#' Havens, T.C. and Bezdek, J.C. (2012): An Efficient Formulation of the
#' Improved Visual Assessment of Cluster Tendency (iVAT) Algorithm, _IEEE
#' Transactions on Knowledge and Data Engineering,_ **24**(5), 813--822.
#' \doi{10.1109/TKDE.2011.33}
#'
#' @source The data was obtained from
#' \url{http://glaros.dtc.umn.edu/gkhome/cluto/cluto/download}
#' @keywords datasets
#' @examples
#' data(Chameleon)
#'
#' plot(chameleon_ds4, cex = .1)
NULL


