library(seriation)
set.seed(0)

context("ser_permutation_vector")

p <- sample(1:10)
sp <- ser_permutation_vector(p, method="valid")

expect_identical(length(sp), 10L)
expect_identical(get_order(sp), p)
expect_identical(get_order(rev(sp)), rev(p))


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
expect_identical(permute(1:10, ser_permutation(1:10)), 1:10)
expect_identical(permute(LETTERS[1:10], ser_permutation(1:10)), LETTERS[1:10])
expect_identical(permute(1:10, ser_permutation(10:1)), 10:1)
expect_identical(permute(LETTERS[1:10], ser_permutation(10:1)), LETTERS[10:1])

expect_error(permute(1:10, ser_permutation(1:11)))

## matrix
m <- matrix(runif(9), ncol=3)
expect_identical(permute(m, ser_permutation(1:3, 3:1)), m[,3:1])
expect_identical(permute(m, ser_permutation(3:1, 3:1)), m[3:1,3:1])

expect_error(permute(m, ser_permutation(1:10, 1:9)))
expect_error(permute(m, ser_permutation(1:9, 1:11)))

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
dend <- as.dendrogram(hclust(d))
expect_equal(dend, permute(dend, order.dendrogram(dend)))
expect_equal(rev(dend), permute(dend, rev(order.dendrogram(dend))))

# chances are that a random order will not be perfect
o <- sample(5)
expect_warning(permute(dend, o))

## hclust
hc <- hclust(d)
expect_equal(hc, permute(hc, get_order(hc)))

### Note: rev for hclust adds labels! (So we only compare merge, height and order)
#expect_equal(rev(hc), permute(hc, rev(get_order(hc))))
expect_equal(rev(hc)[1:3], permute(hc, rev(get_order(hc)))[1:3])

expect_warning(permute(hc, o))

context("permutation_matrix2vector")
pv <- sample(1:100)

## convert into a permutation matrix
pm <- permutation_vector2matrix(pv)

## convert back
expect_identical(permutation_matrix2vector(pm), pv)


