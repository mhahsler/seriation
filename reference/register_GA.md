# Register a Genetic Algorithm Seriation Method

Register a GA-based seriation metaheuristic for use with
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
register_GA()

gaperm_mixedMutation(ismProb = 0.8)
```

## Arguments

- ismProb:

  probability to use
  [`GA::gaperm_ismMutation()`](https://github.com/luca-scr/GA/reference/ga_Mutation.html)
  (inversion) versus
  [`GA::gaperm_simMutation()`](https://github.com/luca-scr/GA/reference/ga_Mutation.html)
  (simple insertion).

## Value

Nothing.

## Details

Registers the method `"GA"` for
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).
This method can be used to optimize any criterion in package seriation.

The GA uses by default the ordered cross-over (OX) operator. For
mutation, the GA uses a mixture of simple insertion and simple inversion
operators. This mixed operator is created using
`seriation::gaperm_mixedMutation(ismProb = .8)`, where `ismProb` is the
probability that the simple insertion mutation operator is used. See
package GA for a description of other available cross-over and mutation
operators for permutations. The appropriate operator functions in GA
start with `gaperm_`.

We warm start the GA using `"suggestions"` given by several heuristics.
Set `"suggestions"` to `NA` to start with a purely random initial
population.

See Example section for available control parameters.

**Note:** Package GA needs to be installed.

## References

Luca Scrucca (2013): GA: A Package for Genetic Algorithms in R. *Journal
of Statistical Software,* **53**(4), 1–37. URL
[doi:10.18637/jss.v053.i04](https://doi.org/10.18637/jss.v053.i04) .

## See also

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`registry_for_seriation_methods`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Author

Michael Hahsler

## Examples

``` r

if (FALSE) { # \dontrun{
register_GA()
get_seriation_method("dist", "GA")

data(SupremeCourt)
d <- as.dist(SupremeCourt)

## optimize for linear seriation criterion (LS)
o <- seriate(d, "GA", criterion = "LS", verbose = TRUE)
pimage(d, o)

## Note that by default the algorithm is already seeded with a LS heuristic.
## This run is no warm start (no suggestions) and increase run to 100
o <- seriate(d, "GA", criterion = "LS", suggestions = NA, run = 100,
  verbose = TRUE)
pimage(d, o)

o <- seriate(d, "GA", criterion = "LS", suggestions = NA, run = 100,
  verbose = TRUE,  )

pimage(d, o)
} # }
```
