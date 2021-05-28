library(seriation)
library(testthat)

context("map")
map <- seriation:::map
v <- 0:10

expect_equal(map(v), seq(0, 1, length.out = length(v)))
expect_equal(map(v, range = c(100,200)), seq(100, 200, length.out = length(v)))
expect_equal(map(v, range = c(200,100)), seq(200, 100, length.out = length(v)))


expect_error(map(v, from.range = c(200,100)))
expect_error(map(v, from.range = c(0, 5, 10)))

expect_equal(map(rep.int(1, 10)), rep(.5, 10))

m <- outer(0:10, 0:10, "+")
expect_equal(map(m), outer(seq(0, 1, length.out = 11), seq(0, 1, length.out = 11), "+") / 2)


context("map_int")
map_int <- seriation:::map_int

expect_identical(map_int(v, range = c(-100, 100)), as.integer(seq(-100, 100, length.out = length(v))))
