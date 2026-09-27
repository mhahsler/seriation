library(seriation)
library(testthat)

map <- seriation:::map
map_int <- seriation:::map_int

test_that("map scales vectors and matrices", {
  v <- 0:10

  expect_equal(map(v), seq(0, 1, length.out = length(v)))
  expect_equal(
    map(v, range = c(100, 200)),
    seq(100, 200, length.out = length(v))
  )
  expect_equal(
    map(v, range = c(200, 100)),
    seq(200, 100, length.out = length(v))
  )
  expect_equal(map(rep.int(1, 10)), rep(0.5, 10))

  m <- outer(0:10, 0:10, "+")
  expected <- outer(
    seq(0, 1, length.out = 11),
    seq(0, 1, length.out = 11),
    "+"
  ) / 2
  expect_equal(map(m), expected)
})

test_that("map validates the source range", {
  v <- 0:10

  expect_error(map(v, from.range = c(200, 100)))
  expect_error(map(v, from.range = c(0, 5, 10)))
})

test_that("map_int returns integer values", {
  v <- 0:10

  expect_identical(
    map_int(v, range = c(-100, 100)),
    as.integer(seq(-100, 100, length.out = length(v)))
  )
})
