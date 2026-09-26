# Neighborhood functions for Seriation Method SA

Definition of different local neighborhood functions for the method
`"SA"` for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
LS_swap(o, pos = sample.int(length(o), 2))

LS_insert(o, pos = sample.int(length(o), 2))

LS_reverse(o, pos = sample.int(length(o), 2))

LS_mixed(o, pos = sample.int(length(o), 2))
```

## Arguments

- o:

  an integer vector with the order

- pos:

  random positions used for the local move.

## Value

returns the new order vector representing the random neighbor.

## Details

Local neighborhood functions are `LS_insert`, `LS_swap`, `LS_reverse`,
and `LS_mix` (1/3 insertion, 1/3 swap and 1/3 reverse). Any neighborhood
function can be defined.

## See also

Other helper:
[`lle()`](http://michael.hahsler.net/seriation/reference/lle.md),
[`uniscale()`](http://michael.hahsler.net/seriation/reference/uniscale.md)
