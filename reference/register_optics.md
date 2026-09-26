# Register Seriation Based on OPTICS

Use ordering points to identify the clustering structure (OPTICS) for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_optics()
```

## Value

Nothing.

## Details

Registers the method `"optics"` for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).
This method applies the OPTICS ordering algorithm implemented in
[`dbscan::optics()`](https://rdrr.io/pkg/dbscan/man/optics.html) to
create an ordering.

**Note:** Package dbscan needs to be installed.

## References

Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Joerg Sander
(1999). OPTICS: Ordering Points To Identify the Clustering Structure.
*ACM SIGMOD international conference on Management of data,* ACM Press,
pp. 49-60.
[doi:10.1145/304181.304187](https://doi.org/10.1145/304181.304187)

## See also

[`dbscan::optics()`](https://rdrr.io/pkg/dbscan/man/optics.html).

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

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
