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

expect_true(round(criterion(d, method="AR_deviations"), 6) - 0.504017 < 1e-10)
## 2.000000 - 1.732051 +  2.236068 - 2.000000 = 0.504017

expect_equal(criterion(d, method="Gradient_raw"), 
	structure(4,names="Gradient_raw"))
## 6 - 2 = 4

expect_true(round(criterion(d, method="Gradient_weighted"), 6) - 3.968119 < 1e-10)
## -1 *(1.000000 - 2.236068 + 1.000000 - 2.000000 + 2.236068 - 2.000000 + 2.000000 - 1.732051 + 1.000000 - 1.732051 + 1.000000 - 2.000000 + 1.732051 - 2.000000 + 2.000000 - 2.236068) 
## = 3.968119

## test stress
expect_equal(round(criterion(d, method="Neumann"), 3), 
	structure(57.275, names="Neumann_stress"))
expect_equal(round(criterion(d, method="Moore"), 3), 
	structure(98.291, names="Moore_stress"))

