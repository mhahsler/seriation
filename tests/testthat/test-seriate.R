### NOTE: disabled snapshot testing since the direction of the order is not defined and randomized
###       for some methods.

library(seriation)
library(testthat)

x <- matrix(
  c(1, 1, 0, 0, 0,
    1, 1, 1, 0, 0,
    0, 0, 1, 1, 1,
    1, 0, 1, 1, 1),
  byrow = TRUE,
  ncol = 5,
  dimnames = list(letters[1:4], LETTERS[1:5])
)

d <- dist(x)

test_that("seriate.dist returns expected results", {
  #local_edition(3) # for snapshot testing

  cat("\n") # for cleaner testthat output
  methods <- list_seriation_methods(kind = "dist")
  os <- sapply(methods, function(m) {
    cat("   -> testing", format(m, width = 13), "... ")
    tm <- system.time(o <- seriate(d, method = m))
    cat("took", formatC(tm[3], digits = 4), "s.\n")
    o
  })
  # make sure they are all the right length
  expect_true(all(sapply(os, length) == nrow(x)))

  # check which methods produce hclusts and which integers
  hclusts <- os[sapply(os, function(x)
    inherits(x, "hclust"))]
  expect_setequal(
    object = names(hclusts),
    expected = c(
      "GW",
      "GW_average",
      "GW_complete",
      "GW_single",
      "GW_ward",
      "HC",
      "HC_average",
      "HC_complete",
      "HC_single",
      "HC_ward",
      "OLO",
      "OLO_average",
      "OLO_complete",
      "OLO_single",
      "OLO_ward"
    )
  )
  integers <- os[sapply(os, is.integer)]
  expect_setequal(
    object = names(integers),
    expected = c(
      "ARSA",
      "BBURCG",
      "BBWRCG",
      "Identity",
      "MDS",
      "MDS_angle",
      "MDS_metric",
      "MDS_nonmetric",
      "QAP_2SUM",
      "QAP_BAR",
      "QAP_Inertia",
      "QAP_LS",
      "R2E",
      "Random",
      "SA",
      "Spectral",
      "Spectral_norm",
      "SPIN_NH",
      "SPIN_STS",
      "TSP",
      "VAT"
    )
  )
  expect_setequal(c(names(hclusts), names(integers)), expected = names(os))

  # check all orders are integers
  ORDERS <-
    sapply(
      X = os,
      FUN = get_order,
      dim = 1,
      simplify = FALSE
    )
  for (n in names(ORDERS))
    expect_type(ORDERS[[!!n]], "integer") # see ?testthat::quasi_label RE "!!"

  # check names - first sanity check before more comprehensive checks
  expect_equal(names(get_order(os$Identity)), letters[1:4])

  # check names of integer seriation vectors correspond to correct integers
  for (n in names(integers)) {
    expect_type(names(ORDERS[[!!n]]), "character") # i.e. not NULL
    expect_mapequal(ORDERS[[!!n]], expected = c(
      a = 1,
      b = 2,
      c = 3,
      d = 4
    ))
  }

  # check $labels of hclust seriation vectors remain in original input order
  for (n in names(hclusts)) {
    expect_equal(hclusts[[!!n]][["labels"]], expected = letters[1:4])
  }

  # check names of get_order(<hclust seriation vector>) equal to ordered labels
  for (n in names(hclusts)) {
    expect_equal(object = names(ORDERS[[!!n]]),
      expected = hclusts[[n]][["labels"]][hclusts[[n]][["order"]]])
  }

  # check str snapshot of some deterministic methods
  deterMethods <- c(
    "BBURCG",
    "BBWRCG",
    "GW",
    "GW_average",
    "GW_complete",
    "GW_single",
    "GW_ward",
    "HC",
    "HC_average",
    "HC_complete",
    "HC_single",
    "HC_ward",
    "Identity",
    "MDS",
    "MDS_angle",
    "MDS_metric",
    "R2E",
    "Spectral",
    "Spectral_norm",
    "VAT"
  )

  #expect_snapshot(str(os[deterMethods]))

})

# check seriate errors for bad dist objects
test_that("negative distances and NAs prompt correct seriate.dist errors", {
  dNeg <- d
  dNeg[1] <- -1
  expect_error(seriate(dNeg), "Negative distances not supported")

  dNA <- d
  dNA[1] <- NA
  expect_error(seriate(dNA), "NAs not allowed in distance matrix x")
})

test_that("dist objects without Diag or Upper attributes can be permuted", {
  # eurodist is an object of class dist from built in R package "datasets"
  expect_s3_class(eurodist, "dist")
  expect_identical(attr(eurodist, "Diag"), NULL)
  expect_identical(attr(eurodist, "Upper"), NULL)
  s <- seriate(eurodist, method = "MDS")
  expect_s3_class(p <- permute(eurodist, order = s), "dist")
  expect_false(attr(p, "Diag")) # permutation adds Diag, is this desirable?
  expect_false(attr(p, "Upper"))
  expect_equal(labels(p), names(get_order(s)))
})

### Stress test to find memory access problems with randomized algorithms
#context("memory stress test")
#replicate(1000, seriate(d, method="bburcg"))
#replicate(1000, seriate(d, method="bbwrcg"))
#replicate(1000, seriate(d, method="arsa"))

test_that("seriate.matrix returns expected results", {
  #local_edition(3) # for snapshot testing

  cat("\n") # for cleaner testthat output
  methods <- list_seriation_methods(kind = "matrix")
  os <- sapply(methods, function(m) {
    cat("   -> testing", format(m, width = 13), "... ")
    tm <- system.time(o <- seriate(x, method = m))
    cat("took", formatC(tm[3], digits = 4), "s.\n")
    o
  }, simplify = FALSE)

  # check number and length of orders
  expect_true(all(sapply(os, length) == 2L))
  expect_true(all(sapply(
    os,
    FUN = function(o2)
      sapply(o2, length)
  ) == c(4L, 5L)))

  x_p <- permute(x, os[[1]]) # BEA method
  expect_equal(x_p, x[get_order(os[[1]], 1), get_order(os[[1]], 2)])

  # check labels
  expect_equal(get_order(os$Identity, 1), c(
    a = 1,
    b = 2,
    c = 3,
    d = 4
  ))
  expect_equal(get_order(os$Identity, 2), c(
    A = 1,
    B = 2,
    C = 3,
    D = 4,
    E = 5
  ))
  expect_equal(get_order(os$Reverse, 2), c(
    A = 5,
    B = 4,
    C = 3,
    D = 2,
    E = 1
  ))

  # check snapshot of some deterministic methods
  #deterMethods <- c("CA", "Identity", "PCA", "PCA_angle", "Reverse")
  #expect_snapshot(str(os[deterMethods]))

})

test_that("seriate.matrix with margin returns expected results", {
  #local_edition(3) # for snapshot testing

  cat("\n") # for cleaner testthat output
  methods <- list_seriation_methods(kind = "matrix")
  os <- sapply(methods, function(m) {
    cat("   -> testing", format(m, width = 13), "... ")
    tm <- system.time(o <- seriate(x, method = m, margin = 2))
    cat("took", formatC(tm[3], digits = 4), "s.\n")
    o
  }, simplify = FALSE)

  expect_true(all(sapply(os, length) == 2L))
  expect_true(all(sapply(
    os,
    FUN = function(o2)
      o2[[1]]
  ) == 1:4))
  expect_true(all(sapply(
    os,
    FUN = function(o2)
      length(o2[[2]]) == 5L
  )))

  x_p <- permute(x, os[[1]], margin = 2)
  expect_equal(x_p, x[, get_order(os[[1]], 2)])
})

test_that("data.frame seriation works as expected", {
  #local_edition(3) # for snapshot testing

  df <- as.data.frame(x)
  o <- seriate(df)
  expect_silent(permute(df, o)) # defaults work with no messages/warnings

  expect_message(
    permute(df, o[1]),
    # DEPRECATED: results in a message
    "permute for data.frames with a single seriation order is now deprecated"
  )

  o <- seriate(df, margin = 1)
  expect_equal(as.integer(o[[2]]), 1:5) # columns left in original order

  oPCA <- seriate(df, method = "PCA")
  #expect_snapshot(permute(df, oPCA))
})
