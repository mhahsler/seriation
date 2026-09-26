# Extracting Order Information from a Permutation Object

Method to get the order information from an object of class
[ser_permutation](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
or
[ser_permutation_vector](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md).
Order information can be extracted as a permutation vector, a vector
containing each object's rank or a permutation matrix.

## Usage

``` r
get_order(x, ...)

# S3 method for class 'ser_permutation_vector'
get_order(x, ...)

# S3 method for class 'ser_permutation'
get_order(x, dim = 1, ...)

# S3 method for class 'hclust'
get_order(x, ...)

# S3 method for class 'dendrogram'
get_order(x, ...)

# S3 method for class 'integer'
get_order(x, ...)

# S3 method for class 'numeric'
get_order(x, ...)

get_rank(x, ...)

get_permutation_matrix(x, ...)
```

## Arguments

- x:

  an object of class
  [ser_permutation](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  or
  [ser_permutation_vector](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md).

- ...:

  further arguments are ignored for `get_order()`. For `get_rank()` and
  for `get_permutation_matrix()` the additional arguments are passed on
  to `get_order()` (e.g., as `dim`).

- dim:

  order information for which dimension should be returned?

## Value

Returns an integer permutation vector/a permutation matrix.

## Details

`get_order()` returns the permutation as an integer vector which
arranges the objects in the seriation order. That is, a vector with the
index of the first, second, \\..., n\\-th object in the order defined by
the permutation. These permutation vectors can directly be used to
reorder objects using subsetting with `"["`. *Note:* In seriation we
usually use these order-based permutation vectors. **Note on names:**
While R's [`order()`](https://rdrr.io/r/base/order.html) returns an
unnamed vector, `get_order()` returns names (if available). The names
are the object label corresponding to the index at that position.
Therefore, the names in the order are in the order after the
permutation.

`get_rank()` returns the seriation as an integer vector containing the
rank/position for each objects after the permutation is applied. That
is, a vector with the position of the first, second, \\..., n\\-th
object after permutation. Note: Use
[`order()`](https://rdrr.io/r/base/order.html) to convert ranks back to
an order.

`get_permutation_matrix()` returns a \\n \times n\\ permutation matrix.

## See also

Other permutation:
[`permutation_vector2matrix()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md),
[`permute()`](http://michael.hahsler.net/seriation/reference/permute.md),
[`reorder.hclust()`](http://michael.hahsler.net/seriation/reference/reorder.hclust.md),
[`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md),
[`ser_permutation()`](http://michael.hahsler.net/seriation/reference/ser_permutation.md),
[`ser_permutation_vector()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)

## Author

Michael Hahsler

## Examples

``` r
## create a random ser_permutation_vector
## Note that ser_permutation_vector is a single permutation vector
x <- structure(1:10, names = paste0("X", 1:10))
o <- sample(x)
o
#>  X4  X2  X7  X9  X1  X6  X3  X5  X8 X10 
#>   4   2   7   9   1   6   3   5   8  10 

p <- ser_permutation_vector(o)
p
#> object of class ‘ser_permutation_vector’, ‘integer’
#> contains a permutation vector of length 10
#> used seriation method: 'unknown'

get_order(p)
#>  X4  X2  X7  X9  X1  X6  X3  X5  X8 X10 
#>   4   2   7   9   1   6   3   5   8  10 
get_rank(p)
#>  X1  X2  X3  X4  X5  X6  X7  X8  X9 X10 
#>   5   2   7   1   8   6   3   9   4  10 
get_permutation_matrix(p)
#>     X4 X2 X7 X9 X1 X6 X3 X5 X8 X10
#> X4   0  0  0  1  0  0  0  0  0   0
#> X2   0  1  0  0  0  0  0  0  0   0
#> X7   0  0  0  0  0  0  1  0  0   0
#> X9   0  0  0  0  0  0  0  0  1   0
#> X1   1  0  0  0  0  0  0  0  0   0
#> X6   0  0  0  0  0  1  0  0  0   0
#> X3   0  0  1  0  0  0  0  0  0   0
#> X5   0  0  0  0  1  0  0  0  0   0
#> X8   0  0  0  0  0  0  0  1  0   0
#> X10  0  0  0  0  0  0  0  0  0   1

## reorder objects using subsetting, the provided permute function or by
## multiplying the with the permutation matrix. We use here
x[get_order(p)]
#>  X4  X2  X7  X9  X1  X6  X3  X5  X8 X10 
#>   4   2   7   9   1   6   3   5   8  10 
permute(x, p)
#>  X4  X2  X7  X9  X1  X6  X3  X5  X8 X10 
#>   4   2   7   9   1   6   3   5   8  10 
drop(get_permutation_matrix(p) %*%  x)
#>  X4  X2  X7  X9  X1  X6  X3  X5  X8 X10 
#>   4   2   7   9   1   6   3   5   8  10 

## ser_permutation contains one permutation vector for each dimension
p2 <- ser_permutation(p, sample(5))
p2
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 2-mode data
#> 
#>   vector length seriation method
#> 1            10          unknown
#> 2             5          unknown

get_order(p2, dim = 2)
#> [1] 3 4 5 1 2
get_rank(p2, dim = 2)
#> [1] 4 5 1 2 3
get_permutation_matrix(p2, dim = 2)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    1    0    0
#> [2,]    0    0    0    1    0
#> [3,]    0    0    0    0    1
#> [4,]    1    0    0    0    0
#> [5,]    0    1    0    0    0
```
