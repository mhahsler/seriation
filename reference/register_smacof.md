# Register Seriation Methods from Package smacof

Registers the `"MDS_smacof"` method for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md)
based on multidimensional scaling using stress majorization and the
corresponding `"smacof_stress0"` criterion implemented in package smacof
(de Leeuw & Mair, 2009).

## Usage

``` r
register_smacof()
```

## Value

Nothing.

## Details

Seriation method `"smacof"` implements stress majorization with several
transformation functions. These functions are passed on as the type
control parameter. We default to `"ratio"`, which together with
`"interval"` performs metric MDS. `"ordinal"` can be used for non-metric
MDS. See
[`smacof::smacofSym()`](https://rdrr.io/pkg/smacof/man/smacofSym.html)
for details on the control parameters.

The corresponding criterion called `"smacof_stress0"` is also
registered. There additional parameter `type` is used to specify the
used transformation function. It should agree with the function used for
seriation. See
[`smacof::stress0()`](https://rdrr.io/pkg/smacof/man/stress0.html) for
details on the stress calculation.

**Note:** Package smacof needs to be installed.

## References

Jan de Leeuw, Patrick Mair (2009). Multidimensional Scaling Using
Majorization: SMACOF in R. *Journal of Statistical Software, 31(3),*
1-30. [doi:10.18637/jss.v031.i03](https://doi.org/10.18637/jss.v031.i03)

## See also

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r
if (FALSE) { # \dontrun{
register_smacof()

get_seriation_method("dist", "MDS_smacof")

d <- dist(random.robinson(20, pre = TRUE))

## use Banded AR form with default clustering (complete-link)
o <- seriate(d, "MDS_smacof", verbose = TRUE)
pimage(d, o)

# recalculate stress for the order
MDS_stress(d, o)

# ordinal MDS. stress needs to be calculated using the correct type with stress0
o <- seriate(d, "MDS_smacof", type = "ordinal", verbose = TRUE)
criterion(d, o, method = "smacof_stress0", type = "ordinal")
} # }
```
