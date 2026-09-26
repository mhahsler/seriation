# Seriation Methods implemented in the R package seriation

This document contains the seriation methods currently implemented in
the R package `seriation`. The methods are organized by methods for
distance matrices and data matrices.

This list was created for the following version of seriation:

``` r

library("seriation")
packageVersion('seriation')
```

    ## [1] '1.5.8.1'

Some additional methods need to be registered.

``` r

register_DendSer()
register_optics()
register_smacof()
register_GA()
register_tsne()
register_umap()
```

Seriation methods on the data `x` can be used in

``` r

seriate(x, method = "Spectral", control = NULL, rep = 1L)
```

The list below shows what seriation criterion (if any) is optimized by
the method.

If the method uses randomization (see randomized below), then `rep` can
be used to run the algorithms several times and return the best
permutation.

Available control parameters can be passed on as a list using the
`control` argument.

### Methods for dissimilarity data (dist)

Seriation methods for dissimilarity data reorder the set of objects in
the data. The methods fall into several groups based on the criteria
they try to optimize, constraints (like dendrograms), and the
algorithmic approach.

#### Dendrogram leaf order

These methods create a dendrogram using hierarchical clustering and then
derive the seriation order from the leaf order in the dendrogram. Leaf
reordering may be applied.

##### DendSer

Dendrogram seriation (Earle and Hurley, 2015).

- optimizes: N/A (specified criterion restricted by dendrogram)
- randomized: FALSE
- registered by: register_DendSer()
- control parameters:

|  | default | help |
|:---|:---|:---|
| h | NULL | an hclust object (optional) |
| method | “complete” | hclust linkage method |
| criterion | NULL | criterion to optimize the dendrogram for |
| DendSer_args | NULL | more arguments are passed on to DendSer (? DendSer) |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### DendSer_ARc

Dendrogram seriation for Anti-Robinson form cost (Earle and Hurley,
2015).

- optimizes:
  [ARc](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#arc)
  (Anti-Robinson form cost restricted by dendrogram)
- randomized: FALSE
- registered by: register_DendSer()
- control parameters:

|  | default | help |
|:---|:---|:---|
| h | NULL | an hclust object (optional) |
| method | “complete” | hclust linkage method |
| criterion | NULL | criterion to optimize the dendrogram for |
| DendSer_args | NULL | more arguments are passed on to DendSer (? DendSer) |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### DendSer_BAR

Dendrogram seriation with BAR (Earle and Hurley, 2015).

- optimizes:
  [BAR](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#bar)
  (banded anti-Robinson form restricted by dendrogram)
- randomized: FALSE
- registered by: register_DendSer()
- control parameters:

|  | default | help |
|:---|:---|:---|
| h | NULL | an hclust object (optional) |
| method | “complete” | hclust linkage method |
| criterion | NULL | criterion to optimize the dendrogram for |
| DendSer_args | NULL | more arguments are passed on to DendSer (? DendSer) |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### DendSer_LPL

Dendrogram seriation for Lazy path length (Earle and Hurley, 2015).

- optimizes:
  [Lazy_path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#lazy_path_length)
  (restricted by dendrogram)
- randomized: FALSE
- registered by: register_DendSer()
- control parameters:

|  | default | help |
|:---|:---|:---|
| h | NULL | an hclust object (optional) |
| method | “complete” | hclust linkage method |
| criterion | NULL | criterion to optimize the dendrogram for |
| DendSer_args | NULL | more arguments are passed on to DendSer (? DendSer) |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### DendSer_PL

Dendrogram seriation for Path length (Earle and Hurley, 2015).

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- registered by: register_DendSer()
- control parameters:

|  | default | help |
|:---|:---|:---|
| h | NULL | an hclust object (optional) |
| method | “complete” | hclust linkage method |
| criterion | NULL | criterion to optimize the dendrogram for |
| DendSer_args | NULL | more arguments are passed on to DendSer (? DendSer) |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### GW

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered by the Gruvaeus and Wainer (1972)
heuristic

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters:

|         | default    | help                                   |
|:--------|:-----------|:---------------------------------------|
| hclust  | NULL       | a precomputed hclust object (optional) |
| linkage | “complete” | hclust method                          |

------------------------------------------------------------------------

##### GW_average

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered by the Gruvaeus and Wainer (1972)
heuristic (avg.link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### GW_complete

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered by the Gruvaeus and Wainer (1972)
heuristic (complete link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### GW_single

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered by the Gruvaeus and Wainer (1972)
heuristic (single link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### GW_ward

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered by the Gruvaeus and Wainer (1972)
heuristic (Ward’s method)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### HC

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering

- optimizes: N/A
- randomized: FALSE
- control parameters:

|         | default    | help                                   |
|:--------|:-----------|:---------------------------------------|
| hclust  | NULL       | a precomputed hclust object (optional) |
| linkage | “complete” | hclust method                          |

------------------------------------------------------------------------

##### HC_average

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering (avg. link).

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### HC_complete

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering (complete link).

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### HC_single

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering (single link)

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### HC_ward

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering (Ward’s method).

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### OLO

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered with optimal leaf ordering
(Bar-Joseph et al., 2001)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters:

|         | default    | help                                   |
|:--------|:-----------|:---------------------------------------|
| hclust  | NULL       | a precomputed hclust object (optional) |
| linkage | “complete” | hclust method                          |

------------------------------------------------------------------------

##### OLO_average

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered with optimal leaf ordering
(Bar-Joseph et al., 2001) (avg. link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### OLO_complete

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered with optimal leaf ordering
(Bar-Joseph et al., 2001) (complete link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### OLO_single

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered with optimal leaf ordering
(Bar-Joseph et al., 2001) (single link)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### OLO_ward

Using the order of the leaf nodes in a dendrogram obtained by
hierarchical clustering and reordered with optimal leaf ordering
(Bar-Joseph et al., 2001) (Ward’s method)

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
  (restricted by dendrogram)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

#### Dimensionality reduction

Find a seriation order by reducing the dimensionality to 1 dimension.
This is typically done by minimizing a stress measure or the
reconstruction error.

##### MDS

Order along the 1D classical metric multidimensional scaling

- optimizes:
  [MDS_stress](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#mds_stress)
  (Euclidean distances)
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| add | FALSE | make the distances Euclidean using an additive constant (see ? cmdscale) |

------------------------------------------------------------------------

##### MDS_angle

Order by the angular order in the 2D MDS projection space split by the
largest gap

- optimizes: N/A
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| add | FALSE | make the distances Euclidean using an additive constant (see ? cmdscale) |

------------------------------------------------------------------------

##### isoMDS

Order along the 1D Kruskal’s non-metric multidimensional scaling

- optimizes:
  [MDS_stress](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#mds_stress)
  (with monotonic transformation)
- randomized: FALSE
- control parameters:

|       | default | help                                                    |
|:------|:--------|:--------------------------------------------------------|
| add   | 1e-09   | small constant to avoid 0 distances                     |
| maxit | 50      | maximum number of iterations                            |
| trace | FALSE   | trace optimization                                      |
| tol   | 0.001   | convergence tolerance                                   |
| p     | 2       | power for Minkowski distance in the configuration space |

------------------------------------------------------------------------

##### isomap

Isometric feature mapping ordination (Tenenbaum, 2000)

- optimizes: N/A (Stress on shortest path distances)
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| k | 30 | number of shortest dissimilarities retained for a point |
| path | “shortest” | method used to estimate the shortest path (“shortest”/“extended”) |

------------------------------------------------------------------------

##### monoMDS

Kruskal’s (1964a,b) non-metric multidimensional scaling (NMDS) using
monotone regression.

- optimizes:
  [MDS_stress](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#mds_stress)
  (Kruskal’s monotone regression stress)
- randomized: TRUE
- control parameters:

|           | default  | help                   |
|:----------|:---------|:-----------------------|
| y         |          | See ? monoMDS for help |
| model     | “global” | N/A                    |
| threshold | 0.8      | N/A                    |
| maxit     | 200      | N/A                    |
| weakties  | TRUE     | N/A                    |
| stress    | 1        | N/A                    |
| scaling   | TRUE     | N/A                    |
| pc        | TRUE     | N/A                    |
| smin      | 1e-04    | N/A                    |
| sfgrmin   | 1e-07    | N/A                    |
| sratmax   | 0.999999 | N/A                    |

------------------------------------------------------------------------

##### metaMDS

Nonmetric Multidimensional Scaling with Stable Solution from Random
Starts.

- optimizes:
  [MDS_stress](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#mds_stress)
  (Kruskal’s monotone regression stress)
- randomized: FALSE
- control parameters:

|               | default   | help                   |
|:--------------|:----------|:-----------------------|
| distance      | “bray”    | see ? metaMDS for help |
| try           | 20        | N/A                    |
| trymax        | 20        | N/A                    |
| engine        | “monoMDS” | N/A                    |
| autotransform | TRUE      | N/A                    |
| noshare       | FALSE     | N/A                    |
| wascores      | TRUE      | N/A                    |
| expand        | TRUE      | N/A                    |
| trace         | 0         | N/A                    |
| plot          | FALSE     | N/A                    |
| previous.best |           | N/A                    |
| verbose       | FALSE     | N/A                    |

------------------------------------------------------------------------

##### Sammon_mapping

Order along the 1D Sammon’s non-linear mapping

- optimizes:
  [MDS_stress](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#mds_stress)
  (scale free, weighted stress called Sammon’s error)
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| add | 1e-09 | small constant to avoid 0 distances |
| niter | 100 | maximum number of iterations |
| trace | FALSE | trace optimization |
| magic | 0.2 | initial value of the step size constant in diagonal Newton method |
| tol | 1e-04 | tolerance for stopping in units of stress |

------------------------------------------------------------------------

##### MDS_smacof

Seriation based on multidimensional scaling using stress majorization
(de Leeuw & Mair, 2009).

- optimizes:
  [smacof_stress0](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#smacof_stress0)
  (MDS stress)
- randomized: FALSE
- registered by: register_smacof()
- control parameters:

|  | default | help |
|:---|:---|:---|
| type | “ratio” | MDS type: “interval”, “ratio”, “ordinal” (nonmetric MDS) |
| init | “torgerson” | start configuration method (“torgerson”/“random”) |
| relax | FALSE | use block relaxation for majorization? |
| modulus | 1 | number of smacof iterations per monotone regression call |
| itmax | 1000 | maximum number of iterations |
| eps | 1e-06 | convergence criterion |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### tsne

Use 1D t-distributed stochastic neighbor embedding (t-SNE) of a distance
matrix to create an order (van der Maaten and Hinton, 2008).

- optimizes: N/A
- randomized: TRUE
- registered by: register_tsne()
- control parameters:

|            | default | help                                                  |
|:-----------|:--------|:------------------------------------------------------|
| max_iter   | 1000    | number of iterations                                  |
| theta      | 0.5     | speed/accuracy trade-off (increase for less accuracy) |
| perplexity | NULL    | perplexity parameter (calculated as n - 1 / 3)        |
| eta        | 100     | learning rate                                         |
| mds        | TRUE    | start from a classical MDS solution                   |
| verbose    | FALSE   | N/A                                                   |

------------------------------------------------------------------------

##### umap

Use 1D Uniform manifold approximation and projection (UMAP) embedding of
the distances to create an order (McInnes and Healy, 2018)

- optimizes: N/A
- randomized: TRUE
- registered by: register_umap()
- control parameters:

|                      | default     | help                      |
|:---------------------|:------------|:--------------------------|
| n_neighbors          | NA          | see ? umap::umap for help |
| n_components         | 1           | N/A                       |
| metric               | “euclidean” | N/A                       |
| n_epochs             | 1000        | N/A                       |
| input                | NA          | N/A                       |
| init                 | “spectral”  | N/A                       |
| min_dist             | 0.1         | N/A                       |
| set_op_mix_ratio     | 1           | N/A                       |
| local_connectivity   | 1           | N/A                       |
| bandwidth            | 1           | N/A                       |
| alpha                | 0.001       | N/A                       |
| gamma                | 1           | N/A                       |
| negative_sample_rate | 5           | N/A                       |
| a                    | NA          | N/A                       |
| b                    | NA          | N/A                       |
| spread               | 1           | N/A                       |
| random_state         | NA          | N/A                       |
| transform_state      | NA          | N/A                       |
| knn                  | NA          | N/A                       |
| knn_repeats          | 1           | N/A                       |
| verbose              | FALSE       | N/A                       |
| umap_learn_args      | NA          | N/A                       |

------------------------------------------------------------------------

#### Optimization

These methods try to optimize a seriation criterion directly, typically
using a heuristic approach.

##### ARSA

Minimize the linear seriation criterion using simulated annealing
(Brusco et al, 2008).

- optimizes:
  [LS](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#ls)
  (Linear seriation criterion)
- randomized: TRUE
- control parameters:

|                   | default | help                                          |
|:------------------|:--------|:----------------------------------------------|
| cool              | 0.5     | cooling factor (smaller means faster cooling) |
| tmin              | 1e-04   | stopping temperature                          |
| swap_to_inversion | 0.5     | probability for swap vs inversion local move  |
| try_multiplier    | 100     | number of local move tries per object         |
| verbose           | FALSE   | N/A                                           |

------------------------------------------------------------------------

##### Enumerate

Enumerate all permutations

- optimizes: N/A (set via control criterion))
- randomized: FALSE
- control parameters:

|           | default             | help                          |
|:----------|:--------------------|:------------------------------|
| criterion | “Gradient_weighted” | Criterion measure to optimize |
| verbose   | FALSE               | N/A                           |

------------------------------------------------------------------------

##### BBURCG

Minimize the unweighted row/column gradient by branch-and-bound (Brusco
and Stahl 2005). This is only feasible for a relatively small number of
objects.

- optimizes:
  [Gradient_raw](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#gradient_raw)
  (Unweighted gradient condition)
- randomized: FALSE
- control parameters:

|         | default | help                                                     |
|:--------|:--------|:---------------------------------------------------------|
| eps     | 0       | Distances need to be at least eps to count as violations |
| verbose | FALSE   | N/A                                                      |

------------------------------------------------------------------------

##### BBWRCG

Minimize the weighted row/column gradient by branch-and-bound (Brusco
and Stahl 2005). This is only feasible for a relatively small number of
objects.

- optimizes:
  [Gradient_weighted](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#gradient_weighted)
  (Weighted gradient condition)
- randomized: FALSE
- control parameters:

|         | default | help                                                     |
|:--------|:--------|:---------------------------------------------------------|
| eps     | 0       | Distances need to be at least eps to count as violations |
| verbose | FALSE   | N/A                                                      |

------------------------------------------------------------------------

##### GA

Use a genetic algorithm to optimize for various criteria.

- optimizes: N/A (specified as parameter criterion)
- randomized: TRUE
- registered by: register_GA()
- control parameters:

|  | default | help |
|:---|:---|:---|
| criterion | “BAR” | criterion to be optimized |
| suggestions | c(“TSP”, “QAP_LS”, “Spectral”) | seed the population with these seriation methods |
| selection | function (object, r = 2/(<object@popSize> \* (<object@popSize> - 1)), q = <2/object@popSize>, …) { | selection operator function |
| crossover | function (object, parents, …) { if (gaControl(“useRcpp”)) gaperm_oxCrossover_Rcpp(objec | crossover operator function |
| mutation | function (object, parent, …) { if (runif(1) \> ismProb) GA::gaperm_simMutation(object, p | mutation operator function |
| pcrossover | 0.8 | probability for crossover |
| pmutation | 0.1 | probability of mutations |
| popSize | 100 | population size |
| maxiter | 1000 | maximum number of generations |
| run | 50 | stop after run generations without improvement |
| parallel | FALSE | use multiple cores? |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### GSA

Minimize a specified seriation measure (criterion) using simulated
annealing. (Hahsler, Hornik, and Buchta, 2023)

- optimizes: N/A (set via control criterion)
- randomized: TRUE
- control parameters:

|  | default | help |
|:---|:---|:---|
| criterion | “Gradient_raw” | Criterion measure to optimize |
| cool | 0.5 | cooling factor (smaller means faster cooling) |
| t_min | 1e-07 | stopping temperature |
| localsearch | “LS_insert” | used local search move function |
| try_multiplier | 5 | number of local move tries per object |
| t0 | NA | initial temperature (if NA then it is estimated) |
| p_initial_accept | 0.01 | Probability to accept a bad move at time 0 (used for t0 estimation) |
| warmstart | “Random” | permutation or seriation method for warmstart |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### SGLS

Improve an existing solution using stochastic greedy local search.
(Hahsler, Hornik, and Buchta, 2023)

- optimizes: N/A (set via control criterion)
- randomized: TRUE
- control parameters:

|             | default        | help                                            |
|:------------|:---------------|:------------------------------------------------|
| criterion   | “Gradient_raw” | Criterion measure to optimize                   |
| init        | “Spectral”     | Start permutation or name of a seriation method |
| max_iter    | NULL           | number of iterations                            |
| localsearch | “LS_insert”    | used local search move function                 |
| verbose     | FALSE          | N/A                                             |

------------------------------------------------------------------------

##### QAP_LS

Quadratic assignment problem formulation for seriation solved using a
simulated annealing solver to minimize the Linear Seriation Problem (LS)
criterion (Hubert and Schultz 1976).

- optimizes:
  [LS](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#ls)
  (Linear seriation criterion)
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

##### QAP_2SUM

Quadratic assignment problem formulation for seriation solved using a
simulated annealing solver to minimize the 2-Sum Problem criterion
(Barnard, Pothen, and Simon 1993).

- optimizes:
  [2SUM](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#sum)
  (2-sum criterion)
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

##### QAP_BAR

Quadratic assignment problem formulation for seriation solved using a
simulated annealing solver to minimize the banded anti-Robinson form
(Hahsler, 2017).

- optimizes:
  [BAR](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#bar)
  (Banded anti-robinson form)
- randomized: TRUE
- control parameters:

|     | default                              | help                       |
|:----|:-------------------------------------|:---------------------------|
| b   | function (n) max(1, floor(n \* 0.2)) | bandwidth (default is 20%) |

------------------------------------------------------------------------

##### QAP_Inertia

Quadratic assignment problem formulation for seriation solved using a
simulated annealing solver to minimize the Inertia criterion (Hahsler,
2017).

- optimizes:
  [Inertia](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#inertia)
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

##### Spectral

Spectral seriation (Ding and He 2004) uses a relaxation to minimize the
2-Sum Problem (Barnard, Pothen, and Simon 1993). It uses the order of
the Fiedler vector of the similarity matrix’s Laplacian.

- optimizes:
  [2SUM](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#sum)
  (2-sum criterion)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### Spectral_norm

Spectral seriation (Ding and He 2004) uses a relaxation to minimize the
2-Sum Problem (Barnard, Pothen, and Simon 1993). It uses the order of
the Fiedler vector of the similarity matrix’s normalized Laplacian.

- optimizes:
  [2SUM](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#sum)
  (2-sum criterion)
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### TSP

Minimize Hamiltonian path length with a TSP solver.

- optimizes:
  [Path_length](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#path_length)
- randomized: TRUE
- control parameters:

|         | default               | help                                 |
|:--------|:----------------------|:-------------------------------------|
| method  | “arbitrary insertion” | used TSP method (see ? solve_TSP)    |
| rep     | 10                    | number of random restarts            |
| two_opt | TRUE                  | use the 2-opt improvement heuristic? |

------------------------------------------------------------------------

#### Other methods

##### Identity

Identity permutation

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### optics

Use ordering points to identify the clustering structure (OPTICS) to
create an order

- optimizes: N/A
- randomized: FALSE
- registered by: register_optics()
- control parameters:

|  | default | help |
|:---|:---|:---|
| eps | Inf | upper limit of the size of the epsilon neighborhood (see ? optics) |
| minPts | 5 | minimum density for dense neighborhoods |

------------------------------------------------------------------------

##### R2E

Rank-two ellipse seriation (Chen 2002)

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### Random

Random permutation

- optimizes: N/A
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

##### Reverse

Reversed identity permutation

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### SPIN_NH

Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005). Neighborhood
algorithm to concentrate low distance values around the diagonal.

- optimizes: N/A (Energy)
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| sigma | c(20, 17, 15, 13, 11, 9, 7, 5, 3, 1) | emphasize local (small alpha) or global (large alpha) structure. |
| step | 5 | iterations to run for each sigma value. |
| W_function | NULL | custom function to create the weight matrix W |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

##### SPIN_STS

Sorting Points Into Neighborhoods (SPIN) (Tsafrir 2005). Side-to-Side
algorithm which tries to push out large distance values.

- optimizes: N/A (Energy)
- randomized: FALSE
- control parameters:

|         | default                         | help                             |
|:--------|:--------------------------------|:---------------------------------|
| step    | 25L                             | iterations to run                |
| nstart  | 10L                             | number of random restarts        |
| X       | function (n) seq(n) - (n + 1)/2 | matrix to calculate the W matrix |
| verbose | FALSE                           | N/A                              |

------------------------------------------------------------------------

##### VAT

Visual assessment of clustering tendency (Bezdek and Hathaway (2002).
Creates an order based on Prim’s algorithm for finding a minimum
spanning tree (MST) in a weighted connected graph representing the
distance matrix. The order is given by the order in which the nodes
(objects) are added to the MST.

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

### Methods for data matrices, tables and data.frames

#### Seriating rows and columns simultaneously

Row and column order influence each other.

##### BEA

Bond Energy Algorithm (BEA; McCormick 1972) to maximize the Measure of
Effectiveness of a non-negative matrix.

- optimizes:
  [ME](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#me)
  (Measure of effectiveness)
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

##### BEA_TSP

Use a TSP to optimize the Measure of Effectiveness (Lenstra 1974).

- optimizes:
  [ME](http://michael.hahsler.net/seriation/articles/seriation_criteria.html#me)
  (Measure of effectiveness)
- randomized: TRUE
- control parameters:

|         | default               | help                                 |
|:--------|:----------------------|:-------------------------------------|
| method  | “arbitrary insertion” | used TSP method (see ? solve_TSP)    |
| rep     | 10                    | number of random restarts            |
| two_opt | TRUE                  | use the 2-opt improvement heuristic? |

------------------------------------------------------------------------

##### BK_unconstrained

Implements the method for binary matrices described in Brower and Kile
(1988). Reorders using the mean row and column position of presences
(1s).

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### CA

This method calculates a correspondence analysis of the matrix and
computes an order according to the scores on a correspondence analysis
dimension (Friendly 2023).

- optimizes: N/A
- randomized: FALSE
- control parameters:

|          | default | help                                          |
|:---------|:--------|:----------------------------------------------|
| dim      | 1L      | CA dimension used for reordering              |
| ca_param | NULL    | List with parameters for the call to ca::ca() |

------------------------------------------------------------------------

#### Seriating rows and columns separately using dissimilarities

##### Heatmap

Calculates distances for rows and columns and then independently applies
the specified seriation method for distances.

- optimizes: N/A
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| dist_fun | list(row = function (x, method = “euclidean”, diag = FALSE, upper = FALSE, p = 2) { if (!is.n | A named list with functions to calculate row and column distances |
| seriation_method | list(row = “OLO_complete”, col = “OLO_complete”) | A named list with row and column seriation methods |
| seriation_control | list(row = NULL, col = NULL) | named list with control parameters for the seriation methods |
| scale | “none” | Scale “rows”, “cols”, or “none” |
| verbose | FALSE | N/A |

------------------------------------------------------------------------

#### Seriate rows using the data matrix

These methods need access to the data matrix instead of dissimilarities
to reorder objects (rows). The same approach can be applied to columns.

##### LLE

Find an order using 1D locally linear embedding.

- optimizes: N/A
- randomized: FALSE
- control parameters:

|     | default | help                              |
|:----|:--------|:----------------------------------|
| k   | 30      | used number of neighbors          |
| reg | 2       | regularization method (see ? lle) |

------------------------------------------------------------------------

##### PCA

Uses the projection of the data on its first principal component to
determine the order.

- optimizes: N/A (Least squares for each dimension (for Euclidean
  distances).)
- randomized: FALSE
- control parameters:

|         | default | help                        |
|:--------|:--------|:----------------------------|
| center  | TRUE    | center the data (mean = 0)? |
| scale   | FALSE   | scale to unit variance?     |
| verbose | FALSE   | FALSE                       |

------------------------------------------------------------------------

##### PCA_angle

Uses the angular order in the 2D PCA projection space split by the
largest gap.

- optimizes: N/A
- randomized: FALSE
- control parameters:

|         | default | help                        |
|:--------|:--------|:----------------------------|
| center  | TRUE    | center the data (mean = 0)? |
| scale   | FALSE   | scale to unit variance?     |
| verbose | FALSE   | FALSE                       |

------------------------------------------------------------------------

##### AOE

Order by the angle of the first two eigenvectors for correlation
matrices (Friendly, 2002)

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### Mean

Reorders rows and columns by row and column means.

- optimizes: N/A
- randomized: FALSE
- control parameters:

|  | default | help |
|:---|:---|:---|
| transformation | NULL | transformation function applied before calculating means (e.g., scale) |

------------------------------------------------------------------------

##### Identity

Identity permutation

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### Reverse

Reversed identity permutation

- optimizes: N/A
- randomized: FALSE
- control parameters: no parameters

------------------------------------------------------------------------

##### Random

Random permutation

- optimizes: N/A
- randomized: TRUE
- control parameters: no parameters

------------------------------------------------------------------------

## References

Arabie, P. and L.J. Hubert (1990): The bond energy algorithm revisited,
*IEEE Transactions on Systems, Man, and Cybernetics,* **20**(1),
268–274. <https://doi.org/10.1109/21.47829>

Bar-Joseph, Z., E. D. Demaine, D. K. Gifford, and T. Jaakkola. (2001):
Fast Optimal Leaf Ordering for Hierarchical Clustering.
*Bioinformatics,* **17**(1), 22–29.
<https://doi.org/10.1093/bioinformatics/17.suppl_1.S22>

Barnard, S. T., A. Pothen, and H. D. Simon (1993): A Spectral Algorithm
for Envelope Reduction of Sparse Matrices. *In Proceedings of the 1993
ACM/IEEE Conference on Supercomputing,* 493–502. Supercomputing ’93. New
York, NY, USA: ACM. <https://ieeexplore.ieee.org/document/1263497>

Bezdek, J.C. and Hathaway, R.J. (2002): VAT: a tool for visual
assessment of (cluster) tendency. *Proceedings of the 2002 International
Joint Conference on Neural Networks (IJCNN ’02),* Volume: 3, 2225–2230.
<https://doi.org/10.1109/IJCNN.2002.1007487>

Brower, J.C. and Kile, K.M. (1988): Sedation of an original data matrix
as applied to paleoecology. *Lethaia,* **21**, 79–93.
<https://doi.org/10.1111/j.1502-3931.1988.tb01756.x>

Brusco, M., Koehn, H.F., and Stahl, S. (2008): Heuristic Implementation
of Dynamic Programming for Matrix Permutation Problems in Combinatorial
Data Analysis. *Psychometrika,* **73**(3), 503–522.
<https://doi.org/10.1007/s11336-007-9049-5>

Brusco, M., and Stahl, S. (2005): *Branch-and-Bound Applications in
Combinatorial Data Analysis.* New York: Springer.
<https://doi.org/10.1007/0-387-28810-4>

Chen, C. H. (2002): Generalized Association Plots: Information
Visualization via Iteratively Generated Correlation Matrices.
*Statistica Sinica,* **12**(1), 7–29.

Ding, C. and Xiaofeng He (2004): Linearized cluster assignment via
spectral ordering. *Proceedings of the Twenty-first International
Conference on Machine learning (ICML ’04)*.
<https://doi.org/10.1145/1015330.1015407>

Climer, S. and Xiongnu Zhang (2006): Rearrangement Clustering: Pitfalls,
Remedies, and Applications, *Journal of Machine Learning Research,*
**7**(Jun), 919–943.

D. Earle, C. B. Hurley (2015): Advances in dendrogram seriation for
application to visualization. *Journal of Computational and Graphical
Statistics,* **24**(1), 1–25.

Friendly, M. (2002): Corrgrams: Exploratory Displays for Correlation
Matrices. *The American Statistician,* **56**(4), 316–324.
<https://doi.org/10.1198/000313002533>

Friendly, M. (2023). *vcdExtra: ‘vcd’ Extensions and Additions*. R
package version 0.8-5, <https://CRAN.R-project.org/package=vcdExtra>.

Gruvaeus, G. and Wainer, H. (1972): Two Additions to Hierarchical
Cluster Analysis, *British Journal of Mathematical and Statistical
Psychology,* **25**, 200–206.
<https://doi.org/10.1111/j.2044-8317.1972.tb00491.x>

Hahsler, M. (2017): An experimental comparison of seriation methods for
one-mode two-way data. *European Journal of Operational Research,*
**257**, 133–143. <https://doi.org/10.1016/j.ejor.2016.08.066>

Hahsler M, Buchta C, Hornik K (2023): *seriation: Infrastructure for
Ordering Objects Using Seriation*. R package version 1.5.0.
<https://doi.org/10.32614/CRAN.package.seriation>

Hahsler, M., Hornik, K., and Buchta, C. (2008): Getting things in order:
An introduction to ithe R package seriation. *Journal of Statistical
Software,* **25**(3), 1–34. <https://doi.org/10.18637/jss.v025.i03>

Hubert, Lawrence, and James Schultz (1976): Quadratic Assignment as a
General Data Analysis Strategy. *British Journal of Mathematical and
Statistical Psychology,* **29**(2). Blackwell Publishing Ltd. 190–241.
<https://doi.org/10.1111/j.2044-8317.1976.tb00714.x>

Hurley, Catherine B. (2004): Clustering Visualizations of
Multidimensional Data. *Journal of Computational and Graphical
Statistics,* **13**(4), 788–806.
<https://doi.org/10.1198/106186004X12425>

Kruskal, J.B. (1964). Nonmetric multidimensional scaling: a numerical
method. *Psychometrika,* **29**, 115–129.

Lenstra, J.K (1974): Clustering a Data Array and the Traveling-Salesman
Problem, *Operations Research,* **22**(2) 413–414.
<https://doi.org/10.1287/opre.22.2.413>

Mair P., De Leeuw J. (2015). Unidimensional scaling. In *Wiley StatsRef:
Statistics Reference Online,* Wiley, New York.
<https://doi.org/10.1002/9781118445112.stat06462.pub2>

McCormick, W.T., P.J. Schweitzer and T.W. White (1972): Problem
decomposition and data reorganization by a clustering technique,
*Operations Research,* **20**(5), 993–1009.
<https://doi.org/10.1287/opre.20.5.993>

McInnes, L and Healy, J, UMAP: Uniform Manifold Approximation and
Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018.

Sammon, J. W. (1969) A non-linear mapping for data structure analysis.
*IEEE Trans. Comput.*, **C-18** 401–409.

Tenenbaum, J.B., de Silva, V. & Langford, J.C. (2000) A global network
framework for nonlinear dimensionality reduction. *Science* **290**,
2319-2323.

Tsafrir, D., Tsafrir, I., Ein-Dor, L., Zuk, O., Notterman, D.A. and
Domany, E. (2005): Sorting points into neighborhoods (SPIN): data
analysis and visualization by ordering distance matrices,
*Bioinformatics,* **21**(10) 2301–8.
<https://doi.org/10.1093/bioinformatics/bti329>

van der Maaten, L.J.P.; Hinton, G.E. (Nov 2008). Visualizing Data Using
t-SNE. *Journal of Machine Learning Research*. **9**: 2579–2605.
