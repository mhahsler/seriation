# Register Seriation Methods from Package DendSer

Register the DendSer dendrogram seriation method and the ARc criterion
(Earle and Hurley, 2015) for use with
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_DendSer()
```

## Value

Nothing.

## Details

Registers the method `"DendSer"` for seriate. DendSer is a fast
heuristic for reordering dendrograms developed by Earle and Hurley
(2015) able to use different criteria.

`control` for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md)
with method `"DendSer"` accepts the following parameters:

- `"h"` or `"method"`: A dendrogram or a method for hierarchical
  clustering (see [hclust](https://rdrr.io/r/stats/hclust.html)).
  Default: complete-link.

- `"criterion"`: A seriation criterion to optimize (see
  `list_criterion_methods("dist")`. Default: `"BAR"` (Banded
  anti-Robinson from with 20% band width).

- `"verbose"`: a logical; print progress information?

- `"DendSer_args"`: additional arguments for
  [`DendSer::DendSer()`](https://rdrr.io/pkg/DendSer/man/DendSer.html).

For convenience, the following methods (for different cost functions)
are also provided:

- `"DendSer_ARc"` (anti-robinson form),

- `"DendSer_BAR"` (banded anti-Robinson form),

- `"DendSer_LPL"` (lazy path length),

- `"DendSer_PL"` (path length).

**Note:** Package DendSer needs to be installed.

## References

D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation for
application to visualization. *Journal of Computational and Graphical
Statistics,* **24**(1), 1–25.

## See also

[`DendSer::DendSer()`](https://rdrr.io/pkg/DendSer/man/DendSer.html)

Other seriation:
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Author

Michael Hahsler based on code by Catherine B. Hurley and Denise Earle

## Examples

``` r

if (FALSE) { # \dontrun{
register_DendSer()
get_seriation_method("dist", "DendSer")

d <- dist(random.robinson(20, pre=TRUE))

## use Banded AR form with default clustering (complete-link)
o <- seriate(d, "DendSer_BAR")
pimage(d, o)

## use different hclust method (Ward) and AR as the cost function for
## dendrogram reordering
o <- seriate(d, "DendSer", control = list(method = "ward.D2", criterion = "AR"))
pimage(d, o)
} # }
```
