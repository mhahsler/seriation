# Class ser_permutation – A Collection of Permutation Vectors for Seriation

The class `ser_permutation` is a collection of permutation vectors (see
class
[ser_permutation_vector](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)),
one for each dimension (mode) of the data to be permuted.

## Usage

``` r
ser_permutation(x, ...)

# S3 method for class 'ser_permutation'
print(x, ...)

# S3 method for class 'ser_permutation'
summary(object, ...)

# S3 method for class 'ser_permutation'
c(..., recursive = FALSE)

# S3 method for class 'ser_permutation'
object[i, ...]
```

## Arguments

- x, object:

  an object of class `ser_permutation_vector` or any object which can be
  converted into a object of class `ser_permutation` (e.g. an integer
  vector).

- ...:

  vectors for further dimensions.

- recursive:

  ignored.

- i:

  index of the dimension(s) to extract.

## Value

An object of class `ser_permutation`.

## See also

Other permutation:
[`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md),
[`permutation_vector2matrix()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md),
[`permute()`](http://michael.hahsler.net/seriation/reference/permute.md),
[`reorder.hclust()`](http://michael.hahsler.net/seriation/reference/reorder.hclust.md),
[`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md),
[`ser_permutation_vector()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)

## Author

Michael Hahsler

## Examples

``` r
o <- ser_permutation(1:5, 10:1)
o
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 2-mode data
#> 
#>   vector length seriation method
#> 1             5          unknown
#> 2            10          unknown

## length (number of dimensions)
length(o)
#> [1] 2

## get permutation vector for 2nd dimension
get_order(o, 2)
#>  [1] 10  9  8  7  6  5  4  3  2  1

## reverse dimensions
o[2:1]
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 2-mode data
#> 
#>   vector length seriation method
#> 1            10          unknown
#> 2             5          unknown

## combine
o <- c(o, ser_permutation(1:15))
o
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 3-mode data
#> 
#>   vector length seriation method
#> 1             5          unknown
#> 2            10          unknown
#> 3            15          unknown

## get an individual permutation
o[[2]]
#> object of class ‘ser_permutation_vector’, ‘integer’
#> contains a permutation vector of length 10
#> used seriation method: 'unknown'

## reverse the order of a permutation
o[[2]] <- rev(o[[2]])
get_order(o,2)
#>  [1]  1  2  3  4  5  6  7  8  9 10
```
