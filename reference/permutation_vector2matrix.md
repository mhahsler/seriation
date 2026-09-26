# Conversion Between Permutation Vector and Permutation Matrix

Converts between permutation vectors and matrices.

## Usage

``` r
permutation_vector2matrix(x)

permutation_matrix2vector(x)
```

## Arguments

- x:

  A permutation vector (any object that can be converted into a
  permutation vector, e.g., a integer vector or a `hclust` object) or a
  matrix representing a permutation. Arguments are checked.

## Value

- `permutation_vector2matrix()`: returns a permutation matrix.

- `permutation_matrix2vector()`: returns the permutation as a integer
  vector.

## See also

Other permutation:
[`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md),
[`permute()`](http://michael.hahsler.net/seriation/reference/permute.md),
[`reorder.hclust()`](http://michael.hahsler.net/seriation/reference/reorder.hclust.md),
[`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md),
[`ser_permutation()`](http://michael.hahsler.net/seriation/reference/ser_permutation.md),
[`ser_permutation_vector()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)

## Author

Michael Hahsler

## Examples

``` r
## create a random permutation vector
pv <- structure(sample(5), names = paste0("X", 1:5))
pv
#> X1 X2 X3 X4 X5 
#>  1  2  3  5  4 

## convert into a permutation matrix
pm <- permutation_vector2matrix(pv)
pm
#>    X1 X2 X3 X4 X5
#> X1  1  0  0  0  0
#> X2  0  1  0  0  0
#> X3  0  0  1  0  0
#> X4  0  0  0  0  1
#> X5  0  0  0  1  0

## convert back
permutation_matrix2vector(pm)
#> X1 X2 X3 X4 X5 
#>  1  2  3  5  4 
```
