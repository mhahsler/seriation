# seriation - Infrastructure for Ordering Objects Using Seriation - R package

[![CRAN version](http://www.r-pkg.org/badges/version/seriation)](https://cran.r-project.org/package=seriation)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/seriation)](https://cran.r-project.org/package=seriation)
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
* __Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/seriation/build/artifacts) or install via `install_github("mhahsler/seriation")` (requires R package `devtools`) 

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
For dissimilarity data:

 *  Branch-and-bound to minimize the unweighted/weighted column gradient 
 *  DendSer - Dendrogram seriation heuristic to optimize various criteria
 *  GA - Genetic algorithm with warm start to optimize various criteria
 *  HC - Hierarchical clustering (single link, avg. link, complete link) 
 *  GW - Hierarchical clustering reordered by Gruvaeus and Wainer heuristic 
 *  OLO - Hierarchical clustering with optimal leaf ordering 
 *  Identity permutation 
 *  MDS - Multidimensional scaling (metric, non-metric, angle) 
 *  SA - Simulated annealing to minimize anti-Robinson events  
 *  TSP - Traveling sales person solver to minimize Hamiltonian path length 
 *  R2E - Rank-two ellipse seriation 
 *  Random permutation
 *  Spectral seriation (unnormalized, normalized) 
 *  SPIN - Sorting points into neighborhoods (neighborhood algorithm, side-to-site algorithm) 
 *  VAT - Visual assessment of clustering tendency ordering 
 *  QAP - Quadratic assignment problem heuristic (2-SUM, linear seriation, inertia, banded anti-Robinson form)
  
For matrices:

 *  BEA - Bond Energy Algorithm to maximize the measure of effectiveness (ME) 
 *  Identity permutation 
 *  PCA - First principal component or angle on the projection on the first two principal components 
 *  Random permutation 
 *  TSP - Traveling sales person solver to maximize ME 

## Further Information

* Development version of [seriation on github](https://github.com/mhahsler/seriation).
* Michael Hahsler, Kurt Hornik and Christian Buchta, [Getting Things in Order: An Introduction to the R Package seriation,](http://dx.doi.org/10.18637/jss.v025.i03) _Journal of Statistical Software,_ 25(3), 2008.
* [Seriation package vignette](http://cran.r-project.org/web/packages/seriation/vignettes/seriation.pdf) with complete examples.
* [Reference manual](http://cran.r-project.org/web/packages/seriation/seriation.pdf)

_Maintainer:_ [Michael Hahsler](http://michael.hahsler.net)
