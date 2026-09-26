# 2D Data Sets used for the CHAMELEON Clustering Algorithm

Several 2D data sets created to evaluate the CHAMELEON clustering
algorithm in the paper by Karypis et al (1999).

## Format

`chameleon_ds4`: The format is a 8,000 x 2 data.frame.

`chameleon_ds5`: The format is a 8,000 x 2 data.frame.

`chameleon_ds7`: The format is a 10,000 x 2 data.frame.

`chameleon_ds8`: The format is a 8,000 x 2 data.frame.

## References

Karypis, G., EH. Han, V. Kumar (1999): CHAMELEON: A Hierarchical
Clustering Algorithm Using Dynamic Modeling, *IEEE Computer,* **32**(8):
68–75. [doi:10.1109/2.781637](https://doi.org/10.1109/2.781637)

## See also

Other data:
[`Irish`](http://michael.hahsler.net/seriation/reference/Irish.md),
[`Munsingen`](http://michael.hahsler.net/seriation/reference/Munsingen.md),
[`Psych24`](http://michael.hahsler.net/seriation/reference/Psych24.md),
[`SupremeCourt`](http://michael.hahsler.net/seriation/reference/SupremeCourt.md),
[`Townships`](http://michael.hahsler.net/seriation/reference/Townships.md),
[`Wood`](http://michael.hahsler.net/seriation/reference/Wood.md),
[`Zoo`](http://michael.hahsler.net/seriation/reference/Zoo.md),
[`create_lines_data()`](http://michael.hahsler.net/seriation/reference/create_lines_data.md),
[`is.robinson()`](http://michael.hahsler.net/seriation/reference/is.robinson.md)

## Examples

``` r
data(Chameleon)

plot(chameleon_ds4, cex = .1)

plot(chameleon_ds5, cex = .1)

plot(chameleon_ds7, cex = .1)

plot(chameleon_ds8, cex = .1)
```
