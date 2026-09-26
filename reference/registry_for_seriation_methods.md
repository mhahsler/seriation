# Registry for Seriation Methods

A registry to manage methods used by
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

## Usage

``` r
registry_seriate

list_seriation_methods(kind, names_only = TRUE)

get_seriation_method(kind, name)

set_seriation_method(
  kind,
  name,
  definition,
  description = NULL,
  control = list(),
  randomized = FALSE,
  optimizes = NA_character_,
  verbose = FALSE,
  ...
)

# S3 method for class 'seriation_method'
print(x, ...)
```

## Arguments

- kind:

  the data type the method works on. For example, `"dist"`, `"matrix"`
  or `"array"`. If missing, then methods for any type are shown.

- names_only:

  logical; return only the method name. `FALSE` returns also the method
  descriptions.

- name:

  the name for the method used to refer to the method in
  [`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

- definition:

  a function containing the method's code.

- description:

  a description of the method. For example, a long name.

- control:

  a list with control arguments and default values.

- randomized:

  logical; does the algorithm use randomization and re-running the
  algorithm several times will lead to different results (see:
  [`seriate_rep()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)).

- optimizes:

  what criterion does the algorithm try to optimize (see:
  [`list_criterion_methods()`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)).

- verbose:

  logical; print a message when a new method is registered.

- ...:

  further information that is stored for the method in the registry.

- x:

  an object of class "seriation_method" to be printed.

## Value

- `list_seriation_method()` result is a vector of character strings with
  the names of the methods. These names are used for methods in
  [`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

- `get_seriation_method()` returns a given method in form of an object
  of class `"seriation_method"`.

## Details

The functions below are convenience function for the registry
`registry_seriate`.

`list_seriation_method()` lists all available methods for a given data
type (`kind`) (e.g., "dist", "matrix"). The result is a vector of
character strings with the method names that can be used in function
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).
If `kind` is missing, then a list of methods is returned.

`get_seriation_method()` returns detailed information for a given method
in form of an object of class `"seriation_method"`. The information
includes a description, parameters and the implementing function.

With `set_seriation_method()` new seriation methods can be added by the
user. The implementing function (`definition`) needs to have the formal
arguments `x, control` and, for arrays and matrices `margin`, where `x`
is the data object and `control` contains a list with additional
information for the method passed on from
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
and `margin` is a vector specifying what dimensions should be seriated.
The implementation has to return a list of objects which can be coerced
into `ser_permutation_vector` objects (e.g., integer vectors). The
elements in the list have to be in corresponding order to the dimensions
of `x`.

## See also

This registry uses
[registry::registry](https://rdrr.io/pkg/registry/man/registry.html).

Other seriation:
[`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md),
[`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md),
[`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md),
[`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md),
[`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md),
[`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md),
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md),
[`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)

## Author

Michael Hahsler

## Examples

``` r
# Registry
registry_seriate
#> An object of class 'registry' with 59 entries.

# List all seriation methods by type
list_seriation_methods()
#> $array
#> [1] "Identity" "Random"   "Reverse" 
#> 
#> $dist
#>  [1] "ARSA"           "BBURCG"         "BBWRCG"         "Enumerate"     
#>  [5] "GSA"            "GW"             "GW_average"     "GW_complete"   
#>  [9] "GW_single"      "GW_ward"        "HC"             "HC_average"    
#> [13] "HC_complete"    "HC_single"      "HC_ward"        "Identity"      
#> [17] "MDS"            "MDS_angle"      "OLO"            "OLO_average"   
#> [21] "OLO_complete"   "OLO_single"     "OLO_ward"       "QAP_2SUM"      
#> [25] "QAP_BAR"        "QAP_Inertia"    "QAP_LS"         "R2E"           
#> [29] "Random"         "Reverse"        "SGD"            "SGLS"          
#> [33] "SPIN_NH"        "SPIN_STS"       "Sammon_mapping" "Spectral"      
#> [37] "Spectral_norm"  "TSP"            "VAT"            "isoMDS"        
#> [41] "isomap"         "metaMDS"        "monoMDS"       
#> 
#> $matrix
#>  [1] "AOE"              "BEA"              "BEA_TSP"          "BK_unconstrained"
#>  [5] "CA"               "Heatmap"          "Identity"         "LLE"             
#>  [9] "Mean"             "PCA"              "PCA_angle"        "Random"          
#> [13] "Reverse"         
#> 

# List methods for matrix seriation
list_seriation_methods("matrix")
#>  [1] "AOE"              "BEA"              "BEA_TSP"          "BK_unconstrained"
#>  [5] "CA"               "Heatmap"          "Identity"         "LLE"             
#>  [9] "Mean"             "PCA"              "PCA_angle"        "Random"          
#> [13] "Reverse"         

get_seriation_method(name = "BEA")
#> name:        BEA
#> kind:        matrix
#> optimizes:   ME (Measure of effectiveness)
#> randomized:  TRUE
#> description: Bond Energy Algorithm (BEA; McCormick 1972) to maximize
#>              the Measure of Effectiveness of a non-negative matrix.
#> control:
#> no parameters
#> 

# Example for defining a new seriation method (reverse identity function for matrix)

# 1. Create the seriation method: Reverse the row order
#    (NA means no seriation is applied to columns)
seriation_method_reverse_rows <- function(x, control = NULL, margin = c(1, 2)) {
    list(rev(seq(nrow(x))), NA)[margin]
}

# 2. Register new method
set_seriation_method("matrix", "Reverse_rows", seriation_method_reverse_rows,
    description = "Reverse identity order", control = list())

list_seriation_methods("matrix")
#>  [1] "AOE"              "BEA"              "BEA_TSP"          "BK_unconstrained"
#>  [5] "CA"               "Heatmap"          "Identity"         "LLE"             
#>  [9] "Mean"             "PCA"              "PCA_angle"        "Random"          
#> [13] "Reverse"          "Reverse_rows"    
get_seriation_method("matrix", "reverse_rows")
#> name:        Reverse_rows
#> kind:        matrix
#> optimizes:   Other
#> randomized:  FALSE
#> description: Reverse identity order
#> control:
#> no parameters
#> 

# 3. Use the new seriation methods
seriate(matrix(1:12, ncol = 3), "reverse_rows")
#> object of class ‘ser_permutation’, ‘list’
#> contains permutation vectors for 2-mode data
#> 
#>   vector length seriation method
#> 1             4     Reverse_rows
#> 2            NA     Reverse_rows
```
