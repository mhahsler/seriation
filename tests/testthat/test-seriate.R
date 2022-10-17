library(seriation)
library(testthat)

x <- matrix(c(
		1,1,0,0,0,
		1,1,1,0,0,
		0,0,1,1,1,
		1,0,1,1,1
		), byrow = TRUE, ncol = 5,
  dimnames = list(1:4, LETTERS[1:5]))

d <- dist(x)

context("seriate_dist")

methods <- list_seriation_methods(kind = "dist")
os <- sapply(methods, function(m) {
  cat("Doing ", m, " ... ")
  tm <- system.time(o <- seriate(d, method = m))
  cat("took ", tm[3],"s.\n")
  o
})

# make sure they are all the right length
expect_true(all(sapply(os, length) == nrow(x)))

# TODO: check labels
#get_order(os$Identity)

# check seriate errors for bad dist objects
test_that("negative distances and NAs prompt correct seriate.dist errors", {
  dNeg <- d
  dNeg[1] <- -1
  expect_error(seriate(dNeg), "Negative distances not supported")

  dNA <- d
  dNA[1] <- NA
  expect_error(seriate(dNA), "NAs not allowed in distance matrix x")
})

### Stress test to find memory access problems with randomized algorithms
#context("memory stress test")
#replicate(1000, seriate(d, method="bburcg"))
#replicate(1000, seriate(d, method="bbwrcg"))
#replicate(1000, seriate(d, method="arsa"))

context("seriate_matrix")

methods <- list_seriation_methods(kind = "matrix")
os <- sapply(methods, function(m) {
  cat("Doing ", m, " ... ")
  tm <- system.time(o <- seriate(x, method = m))
  cat("took ", tm[3],"s.\n")
  o
}, simplify = FALSE)

# check number and length of orders
expect_true(all(sapply(os, length) == 2L))
expect_true(all(sapply(os, FUN = function(o2) sapply(o2, length)) == c(4L, 5L)))

x_p <- permute(x, os[[1]])
expect_equal(x_p, x[get_order(os[[1]], 1), get_order(os[[1]], 2)])

# TODO: check labels
#get_order(os$Identity, 1)
#get_order(os$Identity, 2)
#get_order(os$Reverse, 2)

context("seriate with margin")

methods <- list_seriation_methods(kind = "matrix")
os <- sapply(methods, function(m) {
  cat("Doing ", m, " ... ")
  tm <- system.time(o <- seriate(x, method = m, margin = 2))
  cat("took ", tm[3],"s.\n")
  o
}, simplify = FALSE)

expect_true(all(sapply(os, length) == 2L))
expect_true(all(sapply(os, FUN = function(o2) o2[[1]]) == 1:4))
expect_true(all(sapply(os, FUN = function(o2) length(o2[[2]]) == 5L)))

x_p <- permute(x, os[[1]], margin = 2)
expect_equal(x_p, x[, get_order(os[[1]], 2)])

context("seriate data.frame")
df <- as.data.frame(x)
o <- seriate(df)
permute(df, o)

seriate(df, method = "PCA")

o <- seriate(df, margin = 1)
## DEPRECATED: results in a message
permute(df, o[1])

permute(df, o)
