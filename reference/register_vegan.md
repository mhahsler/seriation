# Register Seriation Methods from Package vegan

Register the `"isomap"`, `"monoMDS"`, and `"metaMDS"` seriation methods
for dissimilarity matrices. These methods use the corresponding
functions from package vegan.

## Usage

``` r
register_vegan()
```

## Value

Nothing.

## Details

**Note:** Package vegan needs to be installed.

## See also

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r
if (FALSE) { # \dontrun{
register_vegan()
list_seriation_methods("dist")
} # }
```
