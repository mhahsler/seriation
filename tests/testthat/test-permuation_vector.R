library(testthat)
library(seriation)
library(dendextend) ## Needed because it redefined all.equal for dendrograms

set.seed(0)

context("ser_permutation_vector")

p <- sample(10)
names(p) <- paste0("X", p)
sp <- ser_permutation_vector(p, method="valid")

expect_identical(length(sp), 10L)
expect_identical(get_order(sp), p)
expect_identical(get_order(rev(sp)), rev(p))
expect_identical(get_rank(sp), structure(order(p), names = names(p)[order(p)]))


expect_error(ser_permutation_vector(c(1:10, 12L), method="invalid"), "Invalid permutation vector!")
expect_error(ser_permutation_vector(c(1:10, 3L), method="invalid"), "Invalid permutation vector!")

context("ser_permutation")

expect_identical(length(ser_permutation(sp)), 1L)
expect_identical(length(ser_permutation(sp, sp)), 2L)

hc <- hclust(dist(runif(10)))
expect_identical(length(ser_permutation(sp, hc)), 2L)
hc <- ser_permutation_vector(hc, method="hc")
expect_identical(length(ser_permutation(sp, hc, sp)), 3L)
expect_identical(length(ser_permutation(ser_permutation(sp), 1:10)), 2L)

context("permute")

## vector
v <- structure(1:10, names = LETTERS[1:10])
expect_identical(permute(v, ser_permutation(1:10)), v[1:10])
expect_identical(permute(LETTERS[1:10], ser_permutation(1:10)), LETTERS[1:10])
expect_identical(permute(v, ser_permutation(10:1)), v[10:1])
expect_identical(permute(LETTERS[1:10], ser_permutation(10:1)), LETTERS[10:1])

expect_error(permute(v, ser_permutation(1:11)))

## matrix
m <- matrix(runif(9), ncol=3, dimnames = list(1:3, LETTERS[1:3]))
expect_identical(permute(m, ser_permutation(1:3, 3:1)), m[,3:1])
expect_identical(permute(m, ser_permutation(3:1, 3:1)), m[3:1,3:1])

expect_error(permute(m, ser_permutation(1:10, 1:9)))
expect_error(permute(m, ser_permutation(1:9, 1:11)))

expect_identical(permute(m, ser_permutation(3:1, 3:1), margin = 1), m[3:1, ])
expect_identical(permute(m, ser_permutation(3:1, 3:1), margin = 2), m[ , 3:1])
expect_identical(permute(m, ser_permutation(3:1), margin = 1), m[3:1, ])
expect_identical(permute(m, ser_permutation(3:1), margin = 2), m[, 3:1])

## data.frame
df <- as.data.frame(m)
expect_identical(permute(df, ser_permutation(1:3, 3:1)), df[,3:1])
expect_identical(permute(df, ser_permutation(3:1, 3:1)), df[3:1,3:1])

## dist
d <- dist(matrix(runif(25), ncol=5))
attr(d, "call") <- NULL   ### permute removes the call attribute
expect_identical(permute(d, ser_permutation(1:5)), d)
### is_equivalent_to ignores attributes
expect_equivalent(permute(d, ser_permutation(5:1)), as.dist(as.matrix(d)[5:1,5:1]))

expect_error(permute(d, ser_permutation(1:8)))

## list
l <- list(a = 1:10, b = letters[1:5], 25)
expect_identical(permute(l, 3:1), rev(l))

## dendrogram
## FIXME: order.dendrogram in stats adds attribute value so I use
## check.attributes = FALSE, but dendrograms use attributes a lot so
## the check may be pointless
dend <- as.dendrogram(hclust(d))

expect_equal(dend, permute(dend, get_order(dend)), ignore_attr = TRUE)
expect_equal(rev(dend), permute(dend, rev(get_order(dend))), ignore_attr = TRUE)

# chances are that a random order will not be perfect
o <- sample(5)
expect_warning(permute(dend, o))

## hclust
hc <- hclust(d)
expect_equal(hc, permute(hc, get_order(hc)))
## Note: rev for hclust adds labels! (So we only compare merge, height and order)
#expect_equal(rev(hc), permute(hc, rev(get_order(hc))))
expect_equal(as.hclust(rev(as.dendrogram(hc)))[1:3],
  permute(hc, rev(get_order(hc)))[1:3])

expect_warning(permute(hc, o))

context("permutation_matrix2vector")
pv <- 1:5
pm <- permutation_vector2matrix(pv)
expect_true(all(diag(pm) == 1))

pv <- sample(1:100)

## convert into a permutation matrix
pm <- permutation_vector2matrix(pv)

## convert back
expect_identical(permutation_matrix2vector(pm), pv)


