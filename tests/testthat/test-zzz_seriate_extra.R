library(seriation)
library(testthat)
### use zzz in the name so it is done as the last test since it
###   registers more methods that should not be tested with the other tests.



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



test_that("t-SNE seriation returns an order", {
  skip_if_not_installed("Rtsne")

  # Note: t-SNE does not work with duplicate entries, which is an issue.
  register_tsne()
  o <- seriate(d, method = "tsne")
  expect_equal(length(o[[1]]), 4L)

  # o <- seriate(x, method = "tsne")
})

test_that("OPTICS seriation returns an order", {
  skip_if_not_installed("dbscan")

  register_optics()
  o <- seriate(d, method = "optics")
  expect_equal(length(o[[1]]), 4L)
})

test_that("GA seriation returns an order", {
  # This is very slow, so only run 10 iterations and skip it on CRAN.
  skip_on_cran()
  skip_if_not_installed("GA")

  register_GA()
  o <- seriate(d, "GA", maxiter = 10, parallel = FALSE, verb = FALSE)
  expect_equal(length(o[[1]]), 4L)
})

test_that("VAE seriation returns orders", {
  # This produces many messages, and Python leaves temporary files that upset
  # CRAN checks. Only run 10 epochs if the test is enabled manually.
  skip("Automatic VAE test is disabled. Run it manually.")
  skip_if_not_installed("keras")

  suppressMessages({
    register_vae()
    o <- seriate(d, "VAE", epochs = 10)
  })

  expect_equal(length(o[[1]]), 4L)

  o <- seriate(x, "VAE", epochs = 10)
  expect_equal(length(o[[1L]]), 4L)
  expect_equal(length(o[[2L]]), 5L)
})
