# Results of 24 Psychological Test for 8th Grade Students

A data set collected by Holzinger and Swineford (1939) which consists of
the results of 24 psychological tests given to 145 seventh and eighth
grade students in a Chicago suburb. This data set contains the
correlation matrix for the 24 test results. The data set was also used
as an example for visualization of cluster analysis by Ling (1973).

## Format

A 24 x 24 correlation matrix.

## References

Holzinger, K. L., Swineford, F. (1939): A study in factor analysis: The
stability of a bi-factor solution. *Supplementary Educational
Monograph,* No. **48**. Chicago: University of Chicago Press.

Ling, R. L. (1973): A computer generated aid for cluster analysis.
*Communications of the ACM,* **16**(6), pp. 355–361.

## See also

Other data:
[`Chameleon`](https://michael.hahsler.net/seriation/reference/Chameleon.md),
[`Irish`](https://michael.hahsler.net/seriation/reference/Irish.md),
[`Munsingen`](https://michael.hahsler.net/seriation/reference/Munsingen.md),
[`SupremeCourt`](https://michael.hahsler.net/seriation/reference/SupremeCourt.md),
[`Townships`](https://michael.hahsler.net/seriation/reference/Townships.md),
[`Wood`](https://michael.hahsler.net/seriation/reference/Wood.md),
[`Zoo`](https://michael.hahsler.net/seriation/reference/Zoo.md),
[`create_lines_data()`](https://michael.hahsler.net/seriation/reference/create_lines_data.md),
[`is.robinson()`](https://michael.hahsler.net/seriation/reference/is.robinson.md)

## Examples

``` r
data("Psych24")

## create a dist object and also get rid of the one negative entry in the
## correlation matrix
d <- as.dist(1 - abs(Psych24))

pimage(d)


## do hclust as in Ling (1973)
hc <- hclust(d, method = "complete")
plot(hc)


pimage(d, hc)


## use seriation
order <- seriate(d, method = "tsp")
#order <- seriate(d, method = "tsp", control = list(method = "concorde"))
pimage(d, order)
```
