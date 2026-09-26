# Class ser_permutation_vector – A Single Permutation Vector for Seriation

The class `ser_permutation_vector` represents a single permutation
vector.

## Usage

``` r
ser_permutation_vector(x, method = NULL)

# S3 method for class 'ser_permutation_vector'
c(..., recursive = FALSE)

# S3 method for class 'ser_permutation_vector'
rev(x)

get_method(x, printable = FALSE)

# S3 method for class 'ser_permutation_vector'
length(x)

# S3 method for class 'ser_permutation_vector'
print(x, ...)

# S3 method for class 'ser_permutation_vector'
summary(object, ...)
```

## Arguments

- x, object:

  an object if class `ser_permutation_vector`. Options for the
  constructor are: (1) an integer permutation vector, (2) an object of
  class [hclust](https://rdrr.io/r/stats/hclust.html), (3) a numeric
  vector with a MDS configuration, or (4) `NA` to indicate a identity
  permutation.

- method:

  a string representing the method used to obtain the permutation
  vector.

- ...:

  further arguments.

- recursive:

  ignored

- printable:

  a logical; prints "unknown" instead of `NULL` for non-existing
  methods.

## Value

The constructor `ser_permutation_vector()` returns an object a
`ser_permutation_vector`

## Details

A permutation vector maps a set of \\n\\ objects \\\\O_1, O_2, ...,
O_n\\\\ onto itself.

**Ordering Representation:** In seriation we represent a permutation
\\\pi\\ as a vector which lists the objects' indices in their permuted
order. This can be seen as replacing the object in position \\i\\ with
the object in position \\\pi(i)\\. For example, the permutation vector
\\\langle3, 1, 2\rangle\\ indicates that in first position is the object
with index 3 then the object with index 1 and finally the object with
index 2. This representation is often called a (re)arrangement or
ordering. The ordering can be extracted from a permutation vector object
via
[`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md).
Such an ordering can be directly used to subset the list of original
objects with `"["` to apply the permutation.

**Rank Representation:** An alternative way to specify a permutation is
via a list of the ranks of the objects after permutation. This
representation is often called a map or substitution. Ranks can be
extracted from a permutation vector using
[`get_rank()`](http://michael.hahsler.net/seriation/reference/get_order.md).

**Permutation Matrix:** Another popular representation is a permutation
matrix which performs permutations using matrix multiplication. A
permutation matrix can be obtained using
[`get_permutation_matrix()`](http://michael.hahsler.net/seriation/reference/get_order.md).

`ser_permutation_vector` objects are usually packed into a
[ser_permutation](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
object which is a collection (a `list`) of \\k\\ permutation vectors for
\\k\\-mode data.

The constructor `ser_permutation_vector()` checks if the permutation
vector is valid (i.e. if all integers occur exactly once).

## See also

Other permutation:
[`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md),
[`permutation_vector2matrix()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md),
[`permute()`](http://michael.hahsler.net/seriation/reference/permute.md),
[`reorder.hclust()`](http://michael.hahsler.net/seriation/reference/reorder.hclust.md),
[`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md),
[`ser_permutation()`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)

## Author

Michael Hahsler

## Examples

``` r
o <- structure(sample(10), names = paste0("X", 1:10))
o
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   1   9   3   2   5   7   6  10   8   4 

p <- ser_permutation_vector(o, "random")
p
#> object of class ‘ser_permutation_vector’, ‘integer’
#> contains a permutation vector of length 10
#> used seriation method: 'random'

## some methods
length(p)
#> [1] 10
get_method(p)
#> [1] "random"
get_order(p)
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   1   9   3   2   5   7   6  10   8   4 
get_rank(p)
#>  X1  X4  X3 X10  X5  X7  X6  X9  X2  X8 
#>   1   4   3  10   5   7   6   9   2   8 
get_permutation_matrix(p)
#>     X1 X2 X3 X4 X5 X6 X7 X8 X9 X10
#> X1   1  0  0  0  0  0  0  0  0   0
#> X2   0  0  0  0  0  0  0  0  1   0
#> X3   0  0  1  0  0  0  0  0  0   0
#> X4   0  1  0  0  0  0  0  0  0   0
#> X5   0  0  0  0  1  0  0  0  0   0
#> X6   0  0  0  0  0  0  1  0  0   0
#> X7   0  0  0  0  0  1  0  0  0   0
#> X8   0  0  0  0  0  0  0  0  0   1
#> X9   0  0  0  0  0  0  0  1  0   0
#> X10  0  0  0  1  0  0  0  0  0   0

r <- rev(p)
r
#> object of class ‘ser_permutation_vector’, ‘integer’
#> contains a permutation vector of length 10
#> used seriation method: 'random'
get_order(r)
#> X10  X9  X8  X7  X6  X5  X4  X3  X2  X1 
#>   4   8  10   6   7   5   2   3   9   1 

## create a symbolic identity permutation vector (with unknown length)
## Note: This can be used to permute an object, but methods
##       like length and get_order are not available.
ip <- ser_permutation_vector(NA)
ip
#> object of class ‘ser_permutation_vector’, ‘integer’
#> contains a permutation vector of length 0
#> used seriation method: 'identity permutation'
```
