library(seriation)


## FIXME add tests for ser_align

set.seed(0)

x <- list(
  a = 1:100,
  b = 100:1,
  c = sample(100),
  d = sample(100)
)


context("ser_dist")

## Default is Spearman
## first two are largest distance (2)
d <- ser_dist(x)
expect_true(all(d >= 0))
expect_equal(d[1], 2)

## x,y interface
d <- ser_dist(x[[1]], x[[2]])
expect_equal(d[1], 2)

## first two are equal with reverse
d <- ser_dist(x, reverse = TRUE)
expect_true(all(d >=0))
expect_equal(d[1], 0)

## Manhattan is 100 times 50 difference
d <- ser_dist(x, method = "Manhattan")
expect_true(all(d >=0))
expect_equal(d[1], 100*50)

d <- ser_dist(x, method = "Manhattan", reverse = TRUE)
expect_true(all(d >=0))
expect_equal(d[1], 0)

## Hamming is 100
d <- ser_dist(x, method = "Hamming")
expect_true(all(d >=0))
expect_equal(d[1], 100)

d <- ser_dist(x, method = "Hamming", reverse = TRUE)
expect_true(all(d >=0))
expect_equal(d[1], 0)

## PPC (reverse has no effect on PPC)
d <- ser_dist(x, method = "PPC")
expect_true(all(d >=0))
expect_equal(d[1], 0)

## test correlations
context("ser_cor")

## Default is Spearman
## sequence with its reverse
co <- ser_cor(x[[1]], x[[2]])
expect_equal(co, rbind(c(1,-1), c(-1,1)))

co <- ser_cor(x)
expect_identical(dim(co), rep(length(x), 2))
expect_true(all(co >=-1 & co <=1))
expect_equivalent(co[1:2,1:2], rbind(c(1,-1), c(-1,1)))

co <- ser_cor(x, reverse = TRUE)
expect_true(all(co >=-1 & co <=1))
expect_equivalent(co[1:2,1:2], rbind(c(1,1), c(1,1)))

### PPC
co <- ser_cor(x, method ="PPC")
expect_true(all(co >=-1 & co <=1))
expect_equivalent(co[1:2,1:2], rbind(c(1,1), c(1,1)))

## test p-value
co <- ser_cor(x, test = TRUE)
expect_equivalent(attr(co, "p-value")[1:2,1:2], matrix(0, nrow=2, ncol=2))

co <- ser_cor(x, reverse = TRUE, test = TRUE)
expect_equivalent(attr(co, "p-value")[1:2,1:2], matrix(0, nrow=2, ncol=2))


