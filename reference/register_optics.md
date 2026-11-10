# Register Seriation Based on OPTICS

Use ordering points to identify the clustering structure (OPTICS) for
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_optics()
```

## Value

Nothing.

## Details

Registers the method `"optics"` for
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md).
This method applies the OPTICS ordering algorithm implemented in
[`dbscan::optics()`](https://michael.hahsler.net/dbscan/reference/optics.html)
to create an ordering.

**Note:** Package dbscan needs to be installed.

## References

Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander
(1999). OPTICS: Ordering Points To Identify the Clustering Structure.
*ACM SIGMOD international conference on Management of data,* ACM Press,
pp. 49-60.
[doi:10.1145/304181.304187](https://doi.org/10.1145/304181.304187)

## See also

[`dbscan::optics()`](https://michael.hahsler.net/dbscan/reference/optics.html).

Other seriation:
[`register_DendSer()`](https://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](https://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_smacof()`](https://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](https://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](https://michael.hahsler.net/seriation/reference/register_umap.md),
[`register_vegan()`](https://michael.hahsler.net/seriation/reference/register_vegan.md),
[`registry_for_seriation_methods`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r

if (FALSE) { # \dontrun{
register_optics()
get_seriation_method("dist", "optics")

d <- dist(random.robinson(50, pre=TRUE, noise=.1))

o <- seriate(d, method = "optics")
pimage(d, o)
} # }
```
