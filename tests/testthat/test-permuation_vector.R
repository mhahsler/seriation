library(testthat)
library(seriation)
suppressMessages(library(dendextend, quietly = True)) ## Needed because it redefined all.equal for dendrograms

set.seed(0)

p <- sample(10)
names(p) <- paste0("X", p)
sp <- ser_permutation_vector(p, method = "valid")

test_that("ser_permutation_vector stores valid orders and ranks", {
  expect_identical(length(sp), 10L)
  expect_identical(get_order(sp), p)
  expect_identical(get_order(rev(sp)), rev(p))
  expect_identical(
    get_rank(sp),
    structure(order(p), names = names(p)[order(p)])
  )

  expect_error(
    ser_permutation_vector(c(1:10, 12L), method = "invalid"),
    "Invalid permutation vector!"
  )
  expect_error(
    ser_permutation_vector(c(1:10, 3L), method = "invalid"),
    "Invalid permutation vector!"
  )
})

test_that("ser_permutation combines permutation representations", {
  expect_identical(length(ser_permutation(sp)), 1L)
  expect_identical(length(ser_permutation(sp, sp)), 2L)

  hc <- hclust(dist(runif(10)))
  expect_identical(length(ser_permutation(sp, hc)), 2L)
  hc <- ser_permutation_vector(hc, method = "hc")
  expect_identical(length(ser_permutation(sp, hc, sp)), 3L)
  expect_identical(
    length(ser_permutation(ser_permutation(sp), 1:10)),
    2L
  )
})

test_that("permute reorders vectors and matrices", {
  v <- structure(1:10, names = LETTERS[1:10])
  expect_identical(permute(v, ser_permutation(1:10)), v[1:10])
  expect_identical(
    permute(LETTERS[1:10], ser_permutation(1:10)),
    LETTERS[1:10]
  )
  expect_identical(permute(v, ser_permutation(10:1)), v[10:1])
  expect_identical(
    permute(LETTERS[1:10], ser_permutation(10:1)),
    LETTERS[10:1]
  )
  expect_error(permute(v, ser_permutation(1:11)))

  m <- matrix(runif(9), ncol = 3, dimnames = list(1:3, LETTERS[1:3]))
  expect_identical(permute(m, ser_permutation(1:3, 3:1)), m[, 3:1])
  expect_identical(permute(m, ser_permutation(3:1, 3:1)), m[3:1, 3:1])
  expect_error(permute(m, ser_permutation(1:10, 1:9)))
  expect_error(permute(m, ser_permutation(1:9, 1:11)))

  expect_identical(
    permute(m, ser_permutation(3:1, 3:1), margin = 1),
    m[3:1, ]
  )
  expect_identical(
    permute(m, ser_permutation(3:1, 3:1), margin = 2),
    m[, 3:1]
  )
  expect_identical(permute(m, ser_permutation(3:1), margin = 1), m[3:1, ])
  expect_identical(permute(m, ser_permutation(3:1), margin = 2), m[, 3:1])
})

test_that("permute reorders data frames, distances, and lists", {
  m <- matrix(runif(9), ncol = 3, dimnames = list(1:3, LETTERS[1:3]))
  df <- as.data.frame(m)
  expect_identical(permute(df, ser_permutation(1:3, 3:1)), df[, 3:1])
  expect_identical(
    permute(df, ser_permutation(3:1, 3:1)),
    df[3:1, 3:1]
  )

  d <- dist(matrix(runif(25), ncol = 5))
  attr(d, "call") <- NULL # permute removes the call attribute
  expect_identical(permute(d, ser_permutation(1:5)), d)
  expect_equal(
    permute(d, ser_permutation(5:1)),
    as.dist(as.matrix(d)[5:1, 5:1]),
    ignore_attr = TRUE
  )
  expect_error(permute(d, ser_permutation(1:8)))

  l <- list(a = 1:10, b = letters[1:5], 25)
  expect_identical(permute(l, 3:1), rev(l))
})

test_that("permute reorders dendrograms and hclust objects", {
  d <- dist(matrix(runif(25), ncol = 5))
  dend <- as.dendrogram(hclust(d))

  # order.dendrogram() adds a value attribute, so attributes are ignored here.
  expect_equal(dend, permute(dend, get_order(dend)), ignore_attr = TRUE)
  expect_equal(
    rev(dend),
    permute(dend, rev(get_order(dend))),
    ignore_attr = TRUE
  )

  # A random order will almost certainly not be perfect.
  o <- sample(5)
  expect_warning(permute(dend, o))

  hc <- hclust(d)
  expect_equal(hc, permute(hc, get_order(hc)))

  # rev() adds labels to hclust, so compare merge, height, and order only.
  expect_equal(
    as.hclust(rev(as.dendrogram(hc)))[1:3],
    permute(hc, rev(get_order(hc)))[1:3]
  )
  expect_warning(permute(hc, o))
})

test_that("permutation matrices convert back to permutation vectors", {
  identity_matrix <- permutation_vector2matrix(1:5)
  expect_true(all(diag(identity_matrix) == 1))

  pv <- sample(1:100)
  pm <- permutation_vector2matrix(pv)
  expect_identical(permutation_matrix2vector(pm), pv)
})
