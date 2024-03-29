library(seriation)

m <- matrix(c(
		1,1,0,0,0,
		1,1,1,0,0,
		0,0,1,1,1,
		1,0,1,1,1
		), byrow=TRUE, ncol=5)

d <- dist(m)
as.matrix(d)

context("criterion")
expect_equal(criterion(d,method="AR_events"), structure(2, names="AR_events"))
## 2

expect_equal(criterion(d,method="Path_length"), structure(4, names="Path_length"))
## 1+2+1=4

expect_equal(criterion(d,method="Lazy_path_length"),
  structure(8, names="Lazy_path_length"))
## (4-1)*1 + (4-2)*2+ (4-3)*1 = 8

expect_true(zapsmall(round(criterion(d, method="AR_deviations"), 6) - 0.504017) == 0)
## 2.000000 - 1.732051 +  2.236068 - 2.000000 = 0.504017

expect_equal(criterion(d, method="Gradient_raw"),
	structure(4,names="Gradient_raw"))
## 6 - 2 = 4

expect_true(zapsmall(round(criterion(d, method="Gradient_weighted"), 6) - 3.968119) == 0)
## -1 *(1.000000 - 2.236068 + 1.000000 - 2.000000 + 2.236068 - 2.000000 + 2.000000 - 1.732051 + 1.000000 - 1.732051 + 1.000000 - 2.000000 + 1.732051 - 2.000000 + 2.000000 - 2.236068)
## = 3.968119

## test stress
expect_equal(round(criterion(d, method="Neumann"), 3),
	structure(7.787, names="Neumann_stress"))
expect_equal(round(criterion(d, method="Moore"), 3),
	structure(11.539, names="Moore_stress"))

expect_equal(criterion(m, method="Neumann"),
	structure(22, names="Neumann_stress"))
expect_equal(criterion(m, method="Moore"),
	structure(44, names="Moore_stress"))


## RGAR
## for w = 2 -> 1/4
## for w = 3 -> 2/8
expect_error(criterion(d, method="RGAR", w=1))
expect_error(criterion(d, method="RGAR", w=4))

expect_equivalent(criterion(d, method="RGAR", pct=0), .25)
expect_equivalent(criterion(d, method="RGAR", w=2), .25)

expect_equivalent(round(criterion(d, method="RGAR", pct=100), 3), .25)
expect_equivalent(round(criterion(d, method="RGAR", w=3), 3), .25)

expect_equivalent(criterion(d, method="RGAR", w=3, relative = FALSE), 2)

### BAR
expect_error(criterion(d, method="BAR", b=0), "Band")
expect_error(criterion(d, method="BAR", b=4), "Band")

# b=1 -> Ham. path length
expect_equivalent(criterion(d, method="BAR", b=1),
  criterion(d, method="Path_length"))
# b = n-1 -> ARc
expect_equivalent(round(criterion(d, method="BAR", b=3), 3), 21.936)

### Cor R
m <- diag(100)

expect_equivalent(criterion(m, method="Cor_R"), 1.0)
expect_equivalent(criterion(m[nrow(m):1,], method="Cor_R"), -1.0)

# this should be close to 0
set.seed(1234)
r <- replicate(100, criterion(m[sample(nrow(m)),], method="Cor_R"))
# hist(r)
expect_true(abs(mean(r)) < 0.1)

# test for data.frame and table

expect_equal(criterion(as.data.frame(m)), criterion(m))
expect_equal(criterion(as.table(m)), criterion(m))
