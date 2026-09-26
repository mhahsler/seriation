# Register Seriation Based on 1D t-SNE

Use t-distributed stochastic neighbor embedding (t-SNE) for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_tsne()
```

## Value

Nothing.

## Details

Registers the method `"tsne"` for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).
This method applies 1D t-SNE to a data matrix or a distance matrix and
extracts the order from the 1D embedding. To speed up the process, an
initial embedding is created using 1D multi-dimensional scaling (MDS) or
principal components analysis (PCA) which is improved by t-SNE.

The `control` parameter `"mds"` or `"pca"` controls if MDS (for
distances) or PCA (for data matrices) is used to create an initial
embedding. See
[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) to learn
about the other available `control` parameters.

Perplexity is automatically set as the minimum between 30 and the number
of observations. It can be also specified using the control parameter
`"perplexity"`.

**Note:** Package Rtsne needs to be installed.

## References

van der Maaten, L.J.P. & Hinton, G.E., 2008. Visualizing
High-Dimensional Data Using t-SNE. *Journal of Machine Learning
Research,* **9**, pp.2579-2605.

## See also

[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r

if (FALSE) { # \dontrun{
register_tsne()

# distances
get_seriation_method("dist", "tsne")

data(SupremeCourt)
d <- as.dist(SupremeCourt)

o <- seriate(d, method = "tsne", verbose = TRUE)
pimage(d, o)

# look at the returned configuration and plot it
attr(o[[1]], "configuration")
plot_config(o)

# the t-SNE results are also available as an attribute (see ? Rtsne::Rtsne)
attr(o[[1]], "model")

## matrix
get_seriation_method("matrix", "tsne")

data("Zoo")
x <- Zoo

x[,"legs"] <- (x[,"legs"] > 0)

# t-SNE does not allow duplicates
x <- x[!duplicated(x), , drop = FALSE]

class <- x$class
label <- rownames(x)
x <- as.matrix(x[,-17])

o <- seriate(x, method = "tsne", eta = 10, verbose = TRUE)
pimage(x, o, prop = FALSE, row_labels = TRUE, col_labels = TRUE)

# look at the row embedding
plot_config(o[[1]], col = class)
} # }
```
