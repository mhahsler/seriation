# Gene Expression Data for Wood Formation in Poplar Trees

A data matrix containing a sample of the normalized gene expression data
for 6 locations in the stem of Popla trees published in the study by
Herzberg et al (2001). The sample of 136 genes selected by Caraux and
Pinloche (2005).

## Format

The format is a 136 x 6 matrix.

## Source

The data was obtained from the Montpellier Bioinformatics Platform

## References

Hertzberg M., H. Aspeborg, J. Schrader, A. Andersson, R.Erlandsson, K.
Blomqvist, R. Bhalerao, M. Uhlen, T. T. Teeri, J. Lundeberg, Bjoern
Sundberg, P. Nilsson and Goeran Sandberg (2001): A transcriptional
roadmap to wood formation, *PNAS,* **98**(25), 14732–14737.

Caraux G. and Pinloche S. (2005): PermutMatrix: a graphical environment
to arrange gene expression profiles in optimal linear order,
*Bioinformatics,* **21**(7) 1280–1281.

## See also

Other data:
[`Chameleon`](http://michael.hahsler.net/seriation/reference/Chameleon.md),
[`Irish`](http://michael.hahsler.net/seriation/reference/Irish.md),
[`Munsingen`](http://michael.hahsler.net/seriation/reference/Munsingen.md),
[`Psych24`](http://michael.hahsler.net/seriation/reference/Psych24.md),
[`SupremeCourt`](http://michael.hahsler.net/seriation/reference/SupremeCourt.md),
[`Townships`](http://michael.hahsler.net/seriation/reference/Townships.md),
[`Zoo`](http://michael.hahsler.net/seriation/reference/Zoo.md),
[`create_lines_data()`](http://michael.hahsler.net/seriation/reference/create_lines_data.md),
[`is.robinson()`](http://michael.hahsler.net/seriation/reference/is.robinson.md)

## Examples

``` r
data(Wood)
head(Wood)
#>                   P          A          B          C          D          E
#> AI161452 -0.7546223 -2.2447910 -2.4157241 -0.8181829  1.0121892  0.8839819
#> AI161500 -2.0621934  0.2127532  0.3556842  0.2219739 -0.6714808  0.3477471
#> AI161513  0.1708342  1.3265617  0.4093247 -1.2003526 -3.3316990 -2.0194944
#> AI161572 -1.1837279 -1.5292043 -2.1512254 -1.0145349  1.1844282 -0.4033869
#> AI161573 -1.8637857 -2.1495779 -2.5108412 -0.8444706  1.4952223 -1.7662259
#> AI161629  1.5917360  1.0212036 -0.1519370 -1.3543136 -2.7099315 -1.3129411
```
