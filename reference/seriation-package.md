# seriation: Infrastructure for Ordering Objects Using Seriation

Infrastructure for ordering objects with an implementation of several
seriation/sequencing/ordination techniques to reorder matrices,
dissimilarity matrices, and dendrograms. Also provides (optimally)
reordered heatmaps, color images and clustering visualizations like
dissimilarity plots, and visual assessment of cluster tendency plots
(VAT and iVAT). Hahsler et al (2008)
[doi:10.18637/jss.v025.i03](https://doi.org/10.18637/jss.v025.i03) .

## Key functions

- Seriation:
  [`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
  [`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md),
  [`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md),
  [`permute()`](http://michael.hahsler.net/seriation/reference/permute.md)

- Visualization:
  [`pimage()`](http://michael.hahsler.net/seriation/reference/pimage.md),
  [`bertinplot()`](http://michael.hahsler.net/seriation/reference/bertinplot.md),
  [`hmap()`](http://michael.hahsler.net/seriation/reference/hmap.md),
  [`dissplot()`](http://michael.hahsler.net/seriation/reference/dissplot.md),
  [`VAT()`](http://michael.hahsler.net/seriation/reference/VAT.md)

## Available seriation methods and criteria

- [A list with the implemented seriation
  methods](https://mhahsler.github.io/seriation/seriation_methods.html)

- [A visual comparison between seriation
  methods](https://mhahsler.github.io/seriation/comparison.html)

- [A list with the implemented seriation
  criteria](https://mhahsler.github.io/seriation/seriation_criteria.html)

## Quickstart guides

- [How to reorder
  heatmaps](https://mhahsler.github.io/seriation/heatmaps.html)

- [How to reorder correlation
  matrices](https://mhahsler.github.io/seriation/correlation_matrix.html)

- [How to evaluate clusters using dissimilarity
  plots](https://mhahsler.github.io/seriation/clustering.html)

## References

Michael Hahsler, Kurt Hornik, and Christian Buchta. Getting things in
order: An introduction to the R package seriation. Journal of
Statistical Software, 25(3):1–34, March 2008.
[doi:10.18637/jss.v025.i03](https://doi.org/10.18637/jss.v025.i03)

## See also

Useful links:

- <https://github.com/mhahsler/seriation>

- <http://michael.hahsler.net/seriation/>

- Report bugs at <https://github.com/mhahsler/seriation/issues>

## Author

**Maintainer**: Michael Hahsler <mhahsler@lyle.smu.edu>
([ORCID](https://orcid.org/0000-0003-2716-1405)) \[copyright holder\]

Authors:

- Michael Hahsler <mhahsler@lyle.smu.edu>
  ([ORCID](https://orcid.org/0000-0003-2716-1405)) \[copyright holder\]

- Christian Buchta \[copyright holder\]

- Kurt Hornik ([ORCID](https://orcid.org/0000-0003-4198-9911))
  \[copyright holder\]

Other contributors:

- David Barnett \[contributor\]

- Michael Brusco \[contributor, copyright holder\]

- Michael Friendly \[contributor\]

- Hans-Friedrich Koehn \[contributor, copyright holder\]

- Fionn Murtagh \[contributor, copyright holder\]

- Stephanie Stahl \[contributor, copyright holder\]
