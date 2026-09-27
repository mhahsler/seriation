library(seriation)
library(testthat)


## FIXME add tests for ser_align

set.seed(0)

x <- list(
  a = 1:100,
  b = 100:1,
  c = sample(100),
  d = sample(100)
)

test_that("ser_dist calculates rank-based dissimilarities", {
  # The default is Spearman. The first two series are equal under reversal.
  d <- ser_dist(x)
  expect_true(all(d >= 0))
  expect_equal(d[1], 0)

  # Without reversal, the first two series have the largest distance.
  d_norev <- ser_dist(x, reverse = FALSE)
  expect_true(all(d_norev >= 0))
  expect_equal(d_norev[1], 2)

  # x, y interface
  d <- ser_dist(x[[1]], x[[2]])
  expect_equal(d[1], 0)
})

test_that("ser_dist supports Manhattan, Hamming, and PPC distances", {
  # Manhattan is 100 times the average difference of 50.
  d <- ser_dist(x, method = "Manhattan", reverse = FALSE)
  expect_true(all(d >= 0))
  expect_equal(d[1], 100 * 50)

  d <- ser_dist(x, method = "Manhattan")
  expect_true(all(d >= 0))
  expect_equal(d[1], 0)

  d <- ser_dist(x, method = "Hamming", reverse = FALSE)
  expect_true(all(d >= 0))
  expect_equal(d[1], 100)

  d <- ser_dist(x, method = "Hamming")
  expect_true(all(d >= 0))
  expect_equal(d[1], 0)

  # Reversal has no effect on PPC.
  d <- ser_dist(x, method = "PPC")
  expect_true(all(d >= 0))
  expect_equal(d[1], 0)
})

test_that("ser_cor calculates correlations with optional reversal", {
  co <- ser_cor(x[[1]], x[[2]], reverse = FALSE)
  expect_equal(co, rbind(c(1, -1), c(-1, 1)))

  co <- ser_cor(x, reverse = FALSE)
  expect_identical(dim(co), rep(length(x), 2))
  expect_true(all(co >= -1 & co <= 1))
  expect_equal(
    co[1:2, 1:2],
    rbind(c(1, -1), c(-1, 1)),
    ignore_attr = TRUE
  )

  co <- ser_cor(x)
  expect_true(all(co >= -1 & co <= 1))
  expect_equal(
    co[1:2, 1:2],
    rbind(c(1, 1), c(1, 1)),
    ignore_attr = TRUE
  )

  co <- ser_cor(x, method = "PPC")
  expect_true(all(co >= -1 & co <= 1))
  expect_equal(
    co[1:2, 1:2],
    rbind(c(1, 1), c(1, 1)),
    ignore_attr = TRUE
  )
})

test_that("ser_cor can calculate p-values", {
  expected <- matrix(0, nrow = 2, ncol = 2)

  co <- ser_cor(x, test = TRUE)
  expect_equal(attr(co, "p-value")[1:2, 1:2], expected, ignore_attr = TRUE)

  co <- ser_cor(x, reverse = TRUE, test = TRUE)
  expect_equal(attr(co, "p-value")[1:2, 1:2], expected, ignore_attr = TRUE)
})
