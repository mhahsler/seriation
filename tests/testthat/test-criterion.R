library(seriation)

m <- matrix(c(
		1,1,0,0,0,
		1,1,1,0,0,
		0,0,1,1,1,
		1,0,1,1,1
		), byrow=TRUE, ncol=5)

d <- dist(m)

test_that("criterion calculates event and path measures", {
  # There are two anti-Robinson events.
  expect_equal(
    criterion(d, method = "AR_events"),
    structure(2, names = "AR_events")
  )

  # Path length: 1 + 2 + 1 = 4.
  expect_equal(
    criterion(d, method = "Path_length"),
    structure(4, names = "Path_length")
  )

  # Lazy path length: (4 - 1) * 1 + (4 - 2) * 2 + (4 - 3) * 1 = 8.
  expect_equal(
    criterion(d, method = "Lazy_path_length"),
    structure(8, names = "Lazy_path_length")
  )

  expect_equal(
    round(criterion(d, method = "AR_deviations"), 6),
    structure(0.504017, names = "AR_deviations")
  )
})

test_that("criterion calculates gradient and stress measures", {
  expect_equal(
    criterion(d, method = "Gradient_raw"),
    structure(4, names = "Gradient_raw")
  )
  expect_equal(
    round(criterion(d, method = "Gradient_weighted"), 6),
    structure(3.968119, names = "Gradient_weighted")
  )

  expect_equal(
    round(criterion(d, method = "Neumann"), 3),
    structure(7.787, names = "Neumann_stress")
  )
  expect_equal(
    round(criterion(d, method = "Moore"), 3),
    structure(11.539, names = "Moore_stress")
  )

  expect_equal(
    criterion(m, method = "Neumann"),
    structure(22, names = "Neumann_stress")
  )
  expect_equal(
    criterion(m, method = "Moore"),
    structure(44, names = "Moore_stress")
  )
})

test_that("RGAR validates and applies its neighborhood width", {
  expect_error(criterion(d, method = "RGAR", w = 1))
  expect_error(criterion(d, method = "RGAR", w = 4))

  # w = 2 gives 1 / 4; w = 3 gives 2 / 8.
  expect_equal(criterion(d, method = "RGAR", pct = 0), 0.25, ignore_attr = TRUE)
  expect_equal(criterion(d, method = "RGAR", w = 2), 0.25, ignore_attr = TRUE)
  expect_equal(
    round(criterion(d, method = "RGAR", pct = 100), 3),
    0.25,
    ignore_attr = TRUE
  )
  expect_equal(
    round(criterion(d, method = "RGAR", w = 3), 3),
    0.25,
    ignore_attr = TRUE
  )
  expect_equal(
    criterion(d, method = "RGAR", w = 3, relative = FALSE),
    2,
    ignore_attr = TRUE
  )
})

test_that("BAR validates and applies its band width", {
  expect_error(criterion(d, method = "BAR", b = 0), "Band")
  expect_error(criterion(d, method = "BAR", b = 4), "Band")

  # b = 1 is Hamiltonian path length; b = n - 1 is ARc.
  expect_equal(
    criterion(d, method = "BAR", b = 1),
    criterion(d, method = "Path_length"),
    ignore_attr = TRUE
  )
  expect_equal(
    round(criterion(d, method = "BAR", b = 3), 3),
    21.936,
    ignore_attr = TRUE
  )
})

test_that("Cor_R measures ordering correlation", {
  m_cor <- diag(100)

  expect_equal(
    criterion(m_cor, method = "Cor_R"),
    1,
    ignore_attr = TRUE
  )
  expect_equal(
    criterion(m_cor[nrow(m_cor):1, ], method = "Cor_R"),
    -1,
    ignore_attr = TRUE
  )

  set.seed(1234)
  random_cor <- replicate(
    100,
    criterion(m_cor[sample(nrow(m_cor)), ], method = "Cor_R")
  )
  expect_lt(abs(mean(random_cor)), 0.1)
})

test_that("criterion supports data frames and tables", {
  expect_equal(criterion(as.data.frame(m)), criterion(m))
  expect_equal(criterion(as.table(m)), criterion(m))
})
