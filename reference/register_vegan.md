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
[`register_DendSer()`](https://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](https://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](https://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](https://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](https://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](https://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)

## Examples

``` r
if (FALSE) { # \dontrun{
register_vegan()
list_seriation_methods("dist")
} # }
```
