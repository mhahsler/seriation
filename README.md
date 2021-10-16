# seriation - Infrastructure for Ordering Objects Using Seriation - R package

[![CRAN version](http://www.r-pkg.org/badges/version/seriation)](https://cran.r-project.org/package=seriation)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/seriation)](https://cran.r-project.org/package=seriation)
[![Travis-CI Build Status](https://travis-ci.org/mhahsler/seriation.svg?branch=master)](https://travis-ci.org/mhahsler/seriation)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mhahsler/seriation?branch=master&svg=true)](https://ci.appveyor.com/project/mhahsler/seriation)

This package provides the infrastructure for ordering objects 
with an implementation of several
[seriation](https://en.wikipedia.org/wiki/Seriation_(archaeology))/sequencing/[ordination](https://en.wikipedia.org/wiki/Ordination_(statistics)) techniques to reorder matrices, dissimilarity
matrices, and dendrograms (see below for a full list). Also provides (optimally) reordered heatmaps, 
color images and clustering visualizations like dissimilarity plots, and
visual assessment of cluster tendency plots (VAT and iVAT).

## Installation

__Stable CRAN version:__ install from within R with
```R
install.packages("seriation")
```
__Current development version:__ Download package from [AppVeyor](https://ci.appveyor.com/project/mhahsler/seriation/build/artifacts) or install from GitHub (needs devtools).
```R 
library("devtools")
install_github("mhahsler/seriation")
```

## Usage

Load library, read data and calculate distances. Then use default seriation.
```R
library(seriation)
data("iris")
x <- as.matrix(iris[-5])
x <- x[sample(1:nrow(x)),]

d <- dist(x)
order <- seriate(d)
order
```

```
object of class ‘ser_permutation’, ‘list’
contains permutation vectors for 1-mode data

  vector length seriation method
1           150             ARSA
```

Compare quality.
```R
rbind(
 random = criterion(d),
 reordered = criterion(d, order)
)
```

```
          AR_events AR_deviations       RGAR Gradient_raw Gradient_weighted Path_length
random       550620    948833.712 0.49938328          741         -1759.954   392.77766
reordered     54846      9426.094 0.04974243       992214       1772123.418    83.95758
            Inertia Least_squares       ME Moore_stress Neumann_stress     2SUM      LS
random    214602194      78852819 291618.0    927570.00     461133.357 29954845 5669489
reordered 356945979      76487641 402332.1     13593.32       5274.093 17810802 4486900
```

## Available Seriation Method

The following methods are available for dissimilarity data:

 *  ARSA - Simulated annealing (linear seriation)   
 *  Branch-and-bound to minimize the unweighted/weighted column gradient 
 *  DendSer - Dendrogram seriation heuristic to optimize various criteria
 *  GA - Genetic algorithm with warm start to optimize various criteria
 *  GW - Hierarchical clustering reordered by Gruvaeus and Wainer heuristic 
 *  HC - Hierarchical clustering (single link, avg. link, complete link) 
 *  Identity permutation 
 *  MDS - Multidimensional scaling (metric, non-metric, angle) 
 *  OLO - Hierarchical clustering with optimal leaf ordering 
 *  OPTICS - Ordering points to identify the clustering structure.
 *  QAP - Quadratic assignment problem heuristic (2-SUM, linear seriation, inertia, banded anti-Robinson form)
 *  R2E - Rank-two ellipse seriation 
 *  Random permutation
 *  Spectral seriation (unnormalized, normalized) 
 *  SPIN - Sorting points into neighborhoods (neighborhood algorithm, side-to-site algorithm)
 *  TSP - Traveling sales person solver to minimize the Hamiltonian path length 
 *  TSNE - Order of the 1D t-distributed stochastic neighbor embedding (t-SNE)
 *  UMAP - Order of the 1D embedding produced by uniform manifold approximation and projection
 *  VAT - Order of the visual assessment of clustering tendency ordering 
  
A detailed comparison of the methods is available in the paper 
[An experimental comparison of seriation methods for one-mode two-way data.](http://dx.doi.org/10.1016/j.ejor.2016.08.066) (read [ preprint](https://michael.hahsler.net/research/misc/EJOR_seriation_2016.pdf)).
  
  
The following methods are available for matrices:

 *  BEA - Bond Energy Algorithm to maximize the measure of effectiveness (ME) 
 *  Identity permutation 
 *  PCA - First principal component or angle on the projection on the first two principal components 
 *  Random permutation 
 *  TSP - Traveling sales person solver to maximize ME 



## References
* [Reference manual for package seriation](https://www.rdocumentation.org/packages/seriation/)
* Michael Hahsler, Kurt Hornik and Christian Buchta, [Getting Things in Order: An Introduction to the R Package seriation,](http://dx.doi.org/10.18637/jss.v025.i03) _Journal of Statistical Software,_ 25(3), 2008.
* Michael Hahsler. [An experimental comparison of seriation methods for one-mode two-way data.](http://dx.doi.org/10.1016/j.ejor.2016.08.066) _European Journal of Operational Research,_ 257:133-143, 2017. (read [preprint](https://michael.hahsler.net/research/misc/EJOR_seriation_2016.pdf))
* [Seriation package vignette](https://cran.r-project.org/package=seriation/vignettes/seriation.pdf) with complete examples.

