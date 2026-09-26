# Registry for Criterion Methods

A registry to manage methods used by
[`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md)
to calculate a criterion value given data and a permutation.

## Usage

``` r
registry_criterion

list_criterion_methods(kind, names_only = TRUE)

get_criterion_method(kind, name)

set_criterion_method(
  kind,
  name,
  fun,
  description = NULL,
  merit = NA,
  control = list(),
  verbose = FALSE,
  ...
)

# S3 method for class 'criterion_method'
print(x, ...)
```

## Arguments

- kind:

  the data type the method works on. For example, `"dist"`, `"matrix"`
  or `"array"`.

- names_only:

  logical; return only the method name. `FALSE` returns also the method
  descriptions.

- name:

  the name for the method used to refer to the method in the function
  [`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md).

- fun:

  a function containing the method's code.

- description:

  a description of the method. For example, a long name.

- merit:

  logical; indicating if the criterion measure is a merit (`TRUE`) or a
  loss (`FALSE`) measure.

- control:

  a list with control arguments and default values.

- verbose:

  logical; print a message when a new method is registered.

- ...:

  further information that is stored for the method in the registry.

- x:

  an object of class "criterion_method" to be printed.

## Value

- `list_criterion_method()` results is a vector of character strings
  with the names of the methods used for
  [`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md).

- `get_criterion_method()` returns a given method in form of an object
  of class `"criterion_method"`.

## Details

All methods below are convenience methods for the registry named
`registry_criterion`.

`list_criterion_method()` lists all available methods for a given data
type (`kind`). The result is a vector of character strings with the
short names of the methods. If `kind` is missing, then a list of methods
is returned.

`get_criterion_method()` returns information (including the implementing
function) about a given method in form of an object of class
`"criterion_method"`.

With `set_criterion_method()` new criterion methods can be added by the
user. The implementing function (`fun`) needs to have the formal
arguments `x, order, ...`, where `x` is the data object, order is an
object of class
[ser_permutation_vector](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
and `...` can contain additional information for the method passed on
from
[`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md).
The implementation has to return the criterion value as a scalar.

## See also

This registry uses
[registry::registry](https://rdrr.io/pkg/registry/man/registry.html).

Other criterion:
[`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md)

## Author

Michael Hahsler

## Examples

``` r
## the registry
registry_criterion
#> An object of class 'registry' with 21 entries.

# List all criterion calculation methods by type
list_criterion_methods()
#> $dist
#>  [1] "2SUM"              "AR_deviations"     "AR_events"        
#>  [4] "BAR"               "Gradient_raw"      "Gradient_weighted"
#>  [7] "Inertia"           "LS"                "Lazy_path_length" 
#> [10] "Least_squares"     "MDS_stress"        "ME"               
#> [13] "Moore_stress"      "Neumann_stress"    "Path_length"      
#> [16] "RGAR"              "Rho"              
#> 
#> $matrix
#> [1] "Cor_R"          "ME"             "Moore_stress"   "Neumann_stress"
#> 

# List methods for matrix
list_criterion_methods("matrix")
#> [1] "Cor_R"          "ME"             "Moore_stress"   "Neumann_stress"

# get more description
list_criterion_methods("matrix", names_only = FALSE)
#> $matrix_Cor_R
#> name:          Cor_R
#> kind:          matrix
#> merit:         TRUE
#> description: Weighted correlation coefficient R: A measure of
#>               effectiveness normalized between -1 and 1 (Deutsch and
#>               Martin, 1971).
#> additional parameters:
#> no parameters
#> 
#> 
#> $matrix_ME
#> name:          ME
#> kind:          matrix
#> merit:         TRUE
#> description: Measure of effectiveness (McCormick, 1972).
#> additional parameters:
#> no parameters
#> 
#> 
#> $matrix_Moore_stress
#> name:          Moore_stress
#> kind:          matrix
#> merit:         FALSE
#> description: Stress criterion (Moore neighborhood) applied to the
#>               reordered matrix (Niermann, 2005).
#> additional parameters:
#> no parameters
#> 
#> 
#> $matrix_Neumann_stress
#> name:          Neumann_stress
#> kind:          matrix
#> merit:         FALSE
#> description: Stress criterion (Neumann neighborhood) applied to the
#>               reordered matrix (Niermann, 2005).
#> additional parameters:
#> no parameters
#> 
#> 

# get a specific method
get_criterion_method(kind = "dist", name = "AR_d")
#> name:          AR_deviations
#> kind:          dist
#> merit:         FALSE
#> description: Anti-Robinson deviations: The number of violations of the
#>               anti-Robinson form weighted by the deviation (Chen,
#>               2002).
#> additional parameters:
#> no parameters
#> 

# Define a new method (sum of the diagonal elements)

## 1. implement a function to calculate the measure
criterion_method_matrix_foo <- function(x, order, ...) {
if(!is.null(order)) x <- permute(x,order)
    sum(diag(x))
}

## 2. Register new method
set_criterion_method("matrix", "DiagSum", criterion_method_matrix_foo,
    description = "Calculated the sum of all diagonal entries", merit = FALSE)

list_criterion_methods("matrix")
#> [1] "Cor_R"          "DiagSum"        "ME"             "Moore_stress"  
#> [5] "Neumann_stress"
get_criterion_method("matrix", "DiagSum")
#> name:          DiagSum
#> kind:          matrix
#> merit:         FALSE
#> description: Calculated the sum of all diagonal entries
#> additional parameters:
#> no parameters
#> 

## 3. use all criterion methods (including the new one)
criterion(matrix(1:9, ncol = 3))
#>          Cor_R        DiagSum             ME   Moore_stress Neumann_stress 
#>    -0.09301487    15.00000000   340.00000000   280.00000000   120.00000000 
```
