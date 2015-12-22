library(seriation)

x <- matrix(c(
		1,1,0,0,0,
		1,1,1,0,0,
		0,0,1,1,1,
		1,0,1,1,1
		), byrow=TRUE, ncol=5) 

d <- dist(x)

context("seriate_dist")

methods <- list_seriation_methods(kind = "dist")
os <- sapply(methods, function(m) {
  cat("Doing ", m, " ... ")
  tm <- system.time(o <- seriate(d, method = m))
  cat("took ", tm[3],"s.\n")
  o
})

### Stress test to find memory access problems with randomized algorithms
#context("memory stress test")
#replicate(1000, seriate(d, method="bburcg"))
#replicate(1000, seriate(d, method="bbwrcg"))
#replicate(1000, seriate(d, method="arsa"))

context("seriate_matrix")

methods <- list_seriation_methods(kind = "matrix")
os <- lapply(methods, function(m) {
  cat("Doing ", m, " ... ")
  tm <- system.time(o <- seriate(x, method = m))
  cat("took ", tm[3],"s.\n")
  o
})
names(os) <- methods

