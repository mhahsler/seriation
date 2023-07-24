### helper to determine the default criterion

get_seriation_kind <- function(x) {
  kind <- class(x)[[1]]
  if (kind %in% c("table", "data.frame"))
    kind <- "matrix"

  kind
}

get_default_criterion <- function(x) {
  kind <- get_seriation_kind(x)

    if (kind == "dist")
      criterion <- "AR_deviations"
    else if (kind == "matrix")
      criterion <- "Moore_stress"
    else stop("Unknown default criterion for type: ", kind)

    criterion
}

get_default_method <- function(x)
  as.list(args(utils::getS3method("seriate", class = class(x)[[1L]])))$method

