# Voting Patterns in the Second Rehnquist U.S. Supreme Court

Contains a (a subset of the) decisions for the stable 8-yr period
1995-2002 of the second Rehnquist Supreme Court. Decisions are
aggregated to the joint probability for disagreement between judges.

## Format

A square, symmetric 9-by-9 matrix with the joint probability for
disagreement.

## References

    Sirovich, L. (2003). A pattern analysis of the second Rehnquist
    U.S. Supreme Court. _Proceedings of the National Academy of Sciences of the United
    States of America,_ **100**, 7432-7437. \doi{10.1073/pnas.1132164100}

## See also

Other data:
[`Chameleon`](http://michael.hahsler.net/seriation/reference/Chameleon.md),
[`Irish`](http://michael.hahsler.net/seriation/reference/Irish.md),
[`Munsingen`](http://michael.hahsler.net/seriation/reference/Munsingen.md),
[`Psych24`](http://michael.hahsler.net/seriation/reference/Psych24.md),
[`Townships`](http://michael.hahsler.net/seriation/reference/Townships.md),
[`Wood`](http://michael.hahsler.net/seriation/reference/Wood.md),
[`Zoo`](http://michael.hahsler.net/seriation/reference/Zoo.md),
[`create_lines_data()`](http://michael.hahsler.net/seriation/reference/create_lines_data.md),
[`is.robinson()`](http://michael.hahsler.net/seriation/reference/is.robinson.md)

## Author

Michael Hahsler

## Examples

``` r
data("SupremeCourt")

# a matrix with joint probability of disagreement
SupremeCourt
#>            Breyer Ginsburg Kennedy OConnor Rehnquist  Scalia  Souter Stevens
#> Breyer    0.00000  0.11966 0.25000 0.20940   0.29915 0.35256 0.11752 0.16239
#> Ginsburg  0.11966  0.00000 0.26790 0.25214   0.30769 0.36966 0.09615 0.14530
#> Kennedy   0.25000  0.26709 0.00000 0.15598   0.12179 0.18803 0.24786 0.32692
#> OConnor   0.20940  0.25214 0.15598 0.00000   0.16239 0.20726 0.22009 0.32906
#> Rehnquist 0.29915  0.30769 0.12179 0.16239   0.00000 0.14316 0.29274 0.40171
#> Scalia    0.35256  0.36966 0.18803 0.20726   0.14316 0.00000 0.33761 0.43803
#> Souter    0.11752  0.09615 0.24790 0.22009   0.29274 0.33761 0.00000 0.16880
#> Stevens   0.16239  0.14530 0.32692 0.32906   0.40171 0.43803 0.16880 0.00000
#> Thomas    0.35897  0.36752 0.17735 0.20513   0.13675 0.06624 0.33120 0.43590
#>            Thomas
#> Breyer    0.35897
#> Ginsburg  0.36752
#> Kennedy   0.17735
#> OConnor   0.20513
#> Rehnquist 0.13675
#> Scalia    0.06624
#> Souter    0.33120
#> Stevens   0.43590
#> Thomas    0.00000

# show judges in original alphabetical order
d <- as.dist(SupremeCourt)
pimage(d, diag = TRUE, upper_tri = TRUE)


# reorder judges using seriation based on similar decisions
o <- seriate(d)
o
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 1-mode data
#> 
#>   vector length seriation method
#> 1             9         Spectral

pimage(d, o, diag = TRUE, upper_tri = TRUE)


# Use optimal leaf ordering (hierarchical clustering with reordering)
# which uses a dendrogram
o <- seriate(d, method = "OLO")
o
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 1-mode data
#> 
#>   vector length seriation method
#> 1             9              OLO

plot(o[[1]])


# Use multi-dimensional scaling and show the configuration
o <- seriate(d, method = "sammon")
o
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 1-mode data
#> 
#>   vector length seriation method
#> 1             9   Sammon_mapping

pimage(d, o, diag = TRUE, upper_tri = TRUE)

plot_config(o[[1]])
```
