# Reorder Dendrograms using Optimal Leaf Ordering

Reorder method for dendrograms for optimal leaf ordering.

## Usage

``` r
# S3 method for class 'hclust'
reorder(x, dist, method = "OLO", ...)
```

## Arguments

- x:

  an object of class `hclust`.

- dist:

  an object of class `dist` with dissimilarities between the objects in
  `x`.

- method:

  a character string with the name of the used measure. Available are:

  - `"OLO"` (optimal leaf ordering; Bar-Joseph et al., 2001) implemented
    in this package and

  - `"GW"` (Gruvaeus and Wainer, 1972) from package gclus.

- ...:

  further arguments are currently ignored.

## Value

A reordered `hclust` object.

## Details

Minimizes the distance between neighboring objects (leaf nodes) in the
dendrogram by flipping the order of subtrees. The algorithm by Gruvaeus
and Wainer is implemented in package gclus (Hurley 2004).

## References

Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola. (2001):
Fast Optimal Leaf Ordering for Hierarchical Clustering.
*Bioinformatics,* **17**(1), 22–29.

Gruvaeus, G. and Wainer, H. (1972): Two Additions to Hierarchical
Cluster Analysis, *British Journal of Mathematical and Statistical
Psychology,* **25**, 200–206.

Hurley, Catherine B. (2004): Clustering Visualizations of
Multidimensional Data. *Journal of Computational and Graphical
Statistics,* **13**(4), 788–806.

## See also

[`gclus::reorder.hclust()`](https://rdrr.io/pkg/gclus/man/hclust.html)

Other permutation:
[`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md),
[`permutation_vector2matrix()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md),
[`permute()`](http://michael.hahsler.net/seriation/reference/permute.md),
[`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md),
[`ser_permutation()`](http://michael.hahsler.net/seriation/reference/ser_permutation.md),
[`ser_permutation_vector()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)

## Author

Michael Hahsler

## Examples

``` r
## cluster European cities by distance
data("eurodist")
d <- as.dist(eurodist)
hc <- hclust(eurodist)

## plot original dendrogram and the reordered dendrograms
plot(hc)

plot(reorder(hc, d, method = "GW"))

plot(reorder(hc, d, method = "OLO"))
```
