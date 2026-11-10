# Register Seriation Based on 1D UMAP

Use uniform manifold approximation and projection (UMAP) to embed the
data on the number line and create a order for
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_umap()
```

## Value

Nothing.

## Details

Registers the method `"umap"` for
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md).
This method applies 1D UMAP to a data matrix or a distance matrix and
extracts the order from the 1D embedding.

Control parameter `n_epochs` can be increased to find a better
embedding.

The returned seriation permutation vector has an attribute named
`embedding` containing the umap embedding.

**Note:** Package umap needs to be installed.

## References

McInnes, L and Healy, J, UMAP: Uniform Manifold Approximation and
Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018.

## See also

[`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html) in umap.

Other seriation:
[`register_DendSer()`](https://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](https://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](https://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](https://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](https://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_vegan()`](https://michael.hahsler.net/seriation/reference/register_vegan.md),
[`registry_for_seriation_methods`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r

if (FALSE) { # \dontrun{
register_umap()

## distances
get_seriation_method("dist", "umap")

data(SupremeCourt)
d <- as.dist(SupremeCourt)

o <- seriate(d, method = "umap", verbose = TRUE)
pimage(d, o)

# look at the returned embedding and plot it
attr(o[[1]], "configuration")
plot_config(o)

## matrix
get_seriation_method("matrix", "umap")

data("Zoo")
Zoo[,"legs"] <- (Zoo[,"legs"] > 0)
x <- as.matrix(Zoo[,-17])
label <- rownames(Zoo)
class <- Zoo$class

o <- seriate(x, method = "umap", verbose = TRUE)
pimage(x, o)

plot_config(o[[1]], col = class)
} # }
```
