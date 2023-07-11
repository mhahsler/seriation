
# <img src="man/figures/logo.svg" align="right" height="139" /> R package seriation - Infrastructure for Ordering Objects Using Seriation

[![CRAN
version](http://www.r-pkg.org/badges/version/seriation)](https://CRAN.R-project.org/package=seriation)
[![stream r-universe
status](https://mhahsler.r-universe.dev/badges/seriation)](https://mhahsler.r-universe.dev/seriation)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/seriation)](https://CRAN.R-project.org/package=seriation)

This package provides the infrastructure for ordering objects with an
implementation of several
[seriation](https://en.wikipedia.org/wiki/Seriation_(archaeology))/sequencing/[ordination](https://en.wikipedia.org/wiki/Ordination_(statistics))
techniques to reorder matrices, dissimilarity matrices, and dendrograms
(see below for a full list). Also provides (optimally) reordered
heatmaps, color images and clustering visualizations like dissimilarity
plots, and visual assessment of cluster tendency plots (VAT and iVAT).

## Available Seriation Method

The following methods are available to reorder dissimilarity data:

**Optimization-based**

- **ARSA** - optimize the linear seriation critreion using simulated
  annealing  
- **Branch-and-bound** to minimize the unweighted/weighted column
  gradient
- **GA** - Genetic algorithm with warm start to optimize any seriation
  criteria
- **GSA** - General simulated annealing to optimize any seriation
  criteria
- **QAP** - Quadratic assignment problem heuristic (2-SUM, linear
  seriation, inertia, banded anti-Robinson form)
- **Spectral** seriation (unnormalized, normalized)
- **TSP** - Traveling sales person solver to minimize the Hamiltonian
  path length

**Multidimensional Scaling**

- **MDS** - classical metric multidimensional scaling
- **Angle in 2D principal coordinates space**
- **isoMDS** - Krusakl’s isoMDS
- **isomap** - Isometric feature mapping ordination
- **monoMDS** - Global and local non-metric multidimensional scaling
  (NMDS).
- **metaMDS** - Nonmetric multidimensional scaling with stable solution
  from random starts
- **Sammon** - Order along the 1D Sammon’s non-linear mapping
- **smacof** - Stress minimization using majorization (metric/nonmetric)
- **TSNE** - Order along the 1D t-distributed stochastic neighbor
  embedding (t-SNE)
- **UMAP** - Order along the 1D embedding produced by uniform manifold
  approximation and projection
- **VAE** - Order along the 1D embedding produced by a variational
  autoencoder.

**Dendrogram-based**

- **DendSer** - Dendrogram seriation heuristic to optimize various
  criteria
- **GW** - Hierarchical clustering reordered by the Gruvaeus and Wainer
  heuristic
- **HC** - Hierarchical clustering (single link, avg. link, complete
  link)
- **OLO** - Hierarchical clustering with optimal leaf ordering

**Other**

- **Identity** permutation
- **OPTICS** - Order of ordering points to identify the clustering
  structure.
- **R2E** - Rank-two ellipse seriation
- **Random** permutation
- **Reverse** order
- **SPIN** - Sorting points into neighborhoods (neighborhood algorithm,
  side-to-site algorithm)
- **VAT** - Order of the visual assessment of clustering tendency

A detailed comparison of the methods is available in the paper [An
experimental comparison of seriation methods for one-mode two-way
data.](http://dx.doi.org/10.1016/j.ejor.2016.08.066) (read the
[preprint](https://michael.hahsler.net/research/paper/EJOR_seriation_2016.pdf)).

The following methods are available to reorder rows and columns of data
matrices:

- **BEA** - Bond Energy Algorithm to maximize the measure of
  effectiveness (ME)
- **BEA_TSP** - TSP to optimize the Measure of Effectiveness
- **CA** - calculates a correspondence analysis of a matrix of
  frequencies (count table) and reorders according to the scores on a
  correspondence analysis dimension
- **Heatmap** - reorders rows and columns independently by calculation
  row/column distances and applying a seriation method (default is
  optimal leaf ordering)
- **Identity** permutation
- **LLE** reorder along a 1D locally linear embedding.
- **Means** - reorders using row and column means.
- **PCA** - orders along the first principal component or angle on the
  projection on the first two principal components
- **Random** permutation
- **Reverse** order
- **TSNE** - Order along the 1D t-distributed stochastic neighbor
  embedding (t-SNE)
- **UMAP** - Order along the 1D embedding produced by uniform manifold
  approximation and projection
- **VAE** - Order along the 1D embedding produced by a variational
  autoencoder.

## Installation

**Stable CRAN version:** Install from within R with

``` r
install.packages("seriation")
```

**Current development version:** Install from
[r-universe.](https://mhahsler.r-universe.dev/seriation)

``` r
install.packages("seriation", repos = "https://mhahsler.r-universe.dev")
```

## Usage

Load library, read data and calculate distances. Then use default
seriation. The used dataset contains the joint probability of
disagreement between Supreme Court Judges during the period from 1995 to
2002.

``` r
library(seriation)
data("SupremeCourt")

SupremeCourt
```

    ##           Breyer Ginsburg Kennedy OConnor Rehnquist Scalia Souter Stevens
    ## Breyer      0.00    0.120    0.25    0.21      0.30  0.353  0.118    0.16
    ## Ginsburg    0.12    0.000    0.27    0.25      0.31  0.370  0.096    0.15
    ## Kennedy     0.25    0.267    0.00    0.16      0.12  0.188  0.248    0.33
    ## OConnor     0.21    0.252    0.16    0.00      0.16  0.207  0.220    0.33
    ## Rehnquist   0.30    0.308    0.12    0.16      0.00  0.143  0.293    0.40
    ## Scalia      0.35    0.370    0.19    0.21      0.14  0.000  0.338    0.44
    ## Souter      0.12    0.096    0.25    0.22      0.29  0.338  0.000    0.17
    ## Stevens     0.16    0.145    0.33    0.33      0.40  0.438  0.169    0.00
    ## Thomas      0.36    0.368    0.18    0.21      0.14  0.066  0.331    0.44
    ##           Thomas
    ## Breyer     0.359
    ## Ginsburg   0.368
    ## Kennedy    0.177
    ## OConnor    0.205
    ## Rehnquist  0.137
    ## Scalia     0.066
    ## Souter     0.331
    ## Stevens    0.436
    ## Thomas     0.000

``` r
d <- as.dist(SupremeCourt)
order <- seriate(d)
order
```

    ## object of class 'ser_permutation', 'list'
    ## contains permutation vectors for 1-mode data
    ## 
    ##   vector length seriation method
    ## 1             9         Spectral

Get the resulting permutation vector.

``` r
get_order(order)
```

    ##    Scalia    Thomas Rehnquist   Kennedy   OConnor    Souter    Breyer  Ginsburg 
    ##         6         9         5         3         4         7         1         2 
    ##   Stevens 
    ##         8

``` r
pimage(d, diag = TRUE, upper = TRUE, main = "judges in original alphabetical order")
pimage(d, order, diag = TRUE, upper = TRUE, main = "judges reordered by seriation")
```

<img src="inst/README_files/seriation-1.png" width="50%" /><img src="inst/README_files/seriation-2.png" width="50%" />

Compare several available quality metrics.

``` r
rbind(alphabetical = criterion(d), seriated = criterion(d, order))
```

    ##              2SUM AR_deviations AR_events BAR Gradient_raw Gradient_weighted
    ## alphabetical  872        10.304        80 1.8            8              0.54
    ## seriated      811         0.064         5 1.1          158             19.76
    ##              Inertia Lazy_path_length Least_squares LS MDS_stress  ME
    ## alphabetical     267              6.9           967 99       0.62  99
    ## seriated         364              4.6           942 86       0.17 101
    ##              Moore_stress Neumann_stress Path_length RGAR   Rho
    ## alphabetical          7.0            3.9         1.8 0.48 0.028
    ## seriated              2.5            1.3         1.1 0.03 0.913

## References

- Michael Hahsler, Kurt Hornik and Christian Buchta, [Getting Things in
  Order: An Introduction to the R Package
  seriation,](http://dx.doi.org/10.18637/jss.v025.i03) *Journal of
  Statistical Software,* 25(3), 2008.
- Michael Hahsler. [An experimental comparison of seriation methods for
  one-mode two-way data.](http://dx.doi.org/10.1016/j.ejor.2016.08.066)
  *European Journal of Operational Research,* 257:133-143, 2017. (read
  the
  [preprint](https://michael.hahsler.net/research/paper/EJOR_seriation_2016.pdf))
- [Reference manual for package
  seriation.](https://www.rdocumentation.org/packages/seriation/)
