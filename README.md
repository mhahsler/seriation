# seriation - Infrastructure for Ordering Objects Using Seriation - R package

[![CRAN version](http://www.r-pkg.org/badges/version/seriation)](http://cran.r-project.org/web/packages/seriation/index.html)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/seriation)](http://cran.r-project.org/web/packages/seriation/index.html)
[![Travis-CI Build Status](https://travis-ci.org/mhahsler/seriation.svg?branch=master)](https://travis-ci.org/mhahsler/seriation)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/seriation?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/seriation)

This package provides the infrastructure for seriation 
with an implementation of several
seriation/sequencing techniques to reorder matrices, dissimilarity
matrices, and dendrograms (see below for a full list). Also provides (optimally) reordered heatmaps, 
color images and clustering visualizations like dissimilarity plots, and
visual assessment of cluster tendency plots (VAT and iVAT).

## Installation

* __Stable CRAN version:__ install from within R.
* __Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/seriation/build/artifacts) or install via `intall_github("mhahsler/seriation")` (requires R package `devtools`) 

## Example

```R
## load library and read data
R> library(seriation)
R> data("iris")
R> x <- as.matrix(iris[-5])
R> x <- x[sample(1:nrow(x)),]

## calculate distances and use default seriation
R> d <- dist(x)
R> order <- seriate(d)
R> order
object of class ‘ser_permutation’, ‘list’
contains permutation vectors for 1-mode data

  vector length seriation method
1           150             ARSA

## compare quality
R> rbind(
+ random = criterion(d),
+ reordered = criterion(d, order)
+ )
          AR_events AR_deviations       RGAR Gradient_raw Gradient_weighted Path_length
random       550620    948833.712 0.49938328          741         -1759.954   392.77766
reordered     54846      9426.094 0.04974243       992214       1772123.418    83.95758
            Inertia Least_squares       ME Moore_stress Neumann_stress     2SUM      LS
random    214602194      78852819 291618.0    927570.00     461133.357 29954845 5669489
reordered 356945979      76487641 402332.1     13593.32       5274.093 17810802 4486900
```

## Available Seriation Methods
For disimilarity data:

 *  2-SUM (QAP) 
 *  Hierarchical clustering (avg. link) 
 *  Hierarchical clustering (avg. link) reordered by Gruvaeus and Wainer heuristic 
 *  Hierarchical clustering (avg. link) with optimal leaf ordering 
 *  Hierarchical clustering (complete link) 
 *  Hierarchical clustering (complete link) reordered by Gruvaeus and Wainer heuristic 
 *  Hierarchical clustering (complete link) with optimal leaf ordering 
 *  Hierarchical clustering reordered by Gruvaeus and Wainer heuristic 
 *  Hierarchical clustering (single link) 
 *  Hierarchical clustering (single link) reordered by Gruvaeus and Wainer heuristic 
 *  Hierarchical clustering (single link) with optimal leaf ordering 
 *  Hierarchical clustering with optimal leaf ordering 
 *  Identity permutation 
 *  Linear Seriation (QAP) 
 *  MDS (angle) 
 *  MDS (metric) 
 *  MDS (non-metric) 
 *  Minimize Anti-Robinson events using simulated annealing 
 *  Minimize Hamiltonian path length with a TSP solver 
 *  Minimize the unweighted row/column gradient by branch-and-bound 
 *  Minimize the weighted row/column gradient by branch-and-bound 
 *  Random permutation 
 *  Rank-two ellipse seriation 
 *  Spectral seriation 
 *  Spectral seriation (normalized) 
 *  SPIN (Neighborhood algorithm) 
 *  SPIN (Side-to-Side algorithm) 
 *  Visual assesment of clustering tendency (VAT) 

For matrices:

 *  Bond Energy Algorithm to maximize ME 
 *  First principal component 
 *  First two principal components (angle) 
 *  Identity permutation 
 *  Random permutation 
 *  TSP to maximize ME 

## Further Information

* Michael Hahsler, Kurt Hornik and Christian Buchta, [Getting Things in Order: An Introduction to the R Package seriation,](http://dx.doi.org/10.18637/jss.v025.i03) _Journal of Statistical Software,_ 25(3), 2008.
* [Seriation package vignette](http://cran.r-project.org/web/packages/seriation/vignettes/seriation.pdf) with complete examples.
* [Reference manual](http://cran.r-project.org/web/packages/seriation/seriation.pdf)


