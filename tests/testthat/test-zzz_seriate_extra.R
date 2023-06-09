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

# this is very slow see we only for 10 iterations
if(seriation:::check_installed("GA", "check")) {
  register_GA()
  o <- seriate(d, "GA", maxiter = 10, parallel = FALSE, verb = F)
  expect_equal(length(o[[1]]), 4L)
}


# Note: tsne does not work with duplicate entries, which is an issue.
if(seriation:::check_installed("Rtsne", "check")) {
  register_tsne()
  o <- seriate(d, method = "tsne")
  expect_equal(length(o[[1]]), 4L)

  #o <- seriate(x, method = "tsne")
}

if(seriation:::check_installed("dbscan", "check")) {
  register_optics()
  o <- seriate(d, method = "optics")
  expect_equal(length(o[[1]]), 4L)
}


# Python (keras) leaves some files in temp and that upsets CRAN
skip_on_cran()

# This produces too many messages
skip()

# only do 10 epochs.
if(seriation:::check_installed("keras", "check")) {
  suppressMessages({
    register_vae()

    o <- seriate(d, "VAE", epochs = 10)
  })

  expect_equal(length(o[[1]]), 4L)

  o <- seriate(x, "VAE", epochs = 10)
  expect_equal(length(o[[1L]]), 4L)
  expect_equal(length(o[[2L]]), 5L)
}
