### NOTE: disabled snapshot testing since the direction of the order is not defined and randomized
###       for some methods.

library(seriation)
library(testthat)

extra_integer <- NULL
extra_hclust <- NULL

if(seriation:::check_installed("DendSer", "check")) {
  register_DendSer()
  extra_hclust <- append(extra_hclust, c("DendSer", "DendSer_ARc",
                                           "DendSer_BAR", "DendSer_LPL",
                                           "DendSer_PL"))
  }

if(seriation:::check_installed("umap", "check")) {
  extra_integer <- append(extra_integer, "umap")
  register_umap()
}

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

test_that("test if seriate.dist returns expected results", {

  cat("\n      dist\n") # for cleaner testthat output
  methods <- list_seriation_methods(kind = "dist")

  ### insufficient data for metaMDS
  methods <- setdiff(methods, "metaMDS")

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
      "OLO_ward",
      extra_hclust
    )
  )
  integers <- os[sapply(os, is.integer)]
  expect_setequal(
    object = names(integers),
    expected = c(
      "ARSA",
      "Enumerate",
      "BBURCG",
      "BBWRCG",
      "Identity",
      "MDS",
      "MDS_angle",
      # "metaMDS",
      "monoMDS",
      "isomap",
      "isoMDS",
      "Sammon_mapping",
      "QAP_2SUM",
      "QAP_BAR",
      "QAP_Inertia",
      "QAP_LS",
      "R2E",
      "Random",
      "Reverse",
      "GSA",
      "SGD",
      "Spectral",
      "Spectral_norm",
      "SPIN_NH",
      "SPIN_STS",
      "TSP",
      "VAT",
      extra_integer
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

  for (o in ORDERS) {
    expect_type(o, "integer")
    expect_mapequal(o, expected = c(
      a = 1,
      b = 2,
      c = 3,
      d = 4
    ))
    expect_type(names(o), "character")
  }

  # check $labels of hclust seriation vectors remain in original input order
  for (n in names(hclusts)) {
    expect_equal(hclusts[[n]][["labels"]], expected = letters[1:4])
  }

  # check names of get_order(<hclust seriation vector>) equal to ordered labels
  for (n in names(hclusts)) {
    expect_equal(object = names(ORDERS[[n]]),
      expected = hclusts[[n]][["labels"]][hclusts[[n]][["order"]]])
  }

  # check snapshot of some deterministic methods
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
    "isoMDS",
    "Sammon_mapping",
    "MDS_angle",
    "R2E",
    "Spectral",
    "Spectral_norm",
    "VAT"
  )

  # recreate with dput(lapply(os[deterMethods], get_order))
  correct <-
    list(
      BBURCG = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      BBWRCG = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      GW = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      GW_average = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      GW_complete = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      GW_single = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      GW_ward = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      HC = structure(1:4, names = c("a",
                                    "b", "c", "d")),
      HC_average = structure(1:4, names = c("a", "b",
                                            "c", "d")),
      HC_complete = structure(1:4, names = c("a", "b",
                                             "c", "d")),
      HC_single = structure(1:4, names = c("a", "b", "c",
                                           "d")),
      HC_ward = structure(1:4, names = c("a", "b", "c", "d")),
      Identity = structure(1:4, names = c("a", "b", "c", "d")),
      MDS = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      isoMDS = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      Sammon_mapping = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      MDS_angle = c(
        a = 1L,
        b = 2L,
        d = 4L,
        c = 3L
      ),
      R2E = c(
        c = 3L,
        d = 4L,
        b = 2L,
        a = 1L
      ),
      Spectral = c(
        c = 3L,
        d = 4L,
        b = 2L,
        a = 1L
      ),
      Spectral_norm = c(c = 3L, d = 4L,
                        b = 2L, a = 1L), VAT = c(c = 3L, d = 4L, b = 2L, a = 1L))

  # some systems may produce the reverse order for some methods!
  for (m in deterMethods)
    expect_true(
      identical(correct[[m]], get_order(os[[m]])) ||
        identical(correct[[m]], rev(get_order(os[[m]]))),
      label = paste("Seriation method", m, "does not return the correct order!\n")
    )
})

# check seriate errors for bad dist objects
test_that("test if negative distances and NAs prompt correct seriate.dist errors", {
  dNeg <- d
  dNeg[1] <- -1
  expect_error(seriate(dNeg), "Negative distances not supported")

  dNA <- d
  dNA[1] <- NA
  expect_error(seriate(dNA), "NAs not allowed in distance matrix x")
})

test_that("test if dist objects without Diag or Upper attributes can be permuted", {
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

test_that("test if seriate.matrix returns expected results", {
  #local_edition(3) # for snapshot testing

  cat("\n      matrix\n") # for cleaner testthat output
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

test_that("test if seriate.matrix with margin returns expected results", {
  #local_edition(3) # for snapshot testing

  cat("\n     matrix with margin\n") # for cleaner testthat output
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

test_that("test if data.frame seriation works as expected", {
  #local_edition(3) # for snapshot testing

  df <- as.data.frame(x)
  o <- seriate(df)
  expect_silent(permute(df, o)) # defaults work with no messages/warnings

  expect_warning(
    permute(df, o[1]),
    # DEPRECATED: results in a message
    "permute for data.frames with a single seriation order is now deprecated"
  )

  o <- seriate(df, margin = 1)
  expect_equal(as.integer(o[[2]]), 1:5) # columns left in original order

  oPCA <- seriate(df, method = "PCA")
  #expect_snapshot(permute(df, oPCA))
})

