# A Comparison of Seriation Methods

## Introduction

This document compares the seriation methods available in the package
`seriation` using a sample of 30 flowers from the popular Iris dataset
and randomizes the order of the objects.

``` r

set.seed(1234)
library("seriation")

data("iris")
x <- as.matrix(iris[sample(nrow(iris), 30), -5])
d <- dist(x)
```

## Distance seriation

We first register more seriation methods. Some of these methods require
installing additional packages.

``` r

register_DendSer()
register_optics()
register_smacof()
```

The following methods will be used (a few slow methods are skipped).

``` r

methods <- sort(list_seriation_methods("dist"))
methods <- setdiff(methods, c("BBURCG", "BBWRCG", "Enumerate", "GSA", "SGD", "SGLS"))
methods 
```

    ##  [1] "ARSA"           "DendSer"        "DendSer_ARc"    "DendSer_BAR"   
    ##  [5] "DendSer_LPL"    "DendSer_PL"     "GW"             "GW_average"    
    ##  [9] "GW_complete"    "GW_single"      "GW_ward"        "HC"            
    ## [13] "HC_average"     "HC_complete"    "HC_single"      "HC_ward"       
    ## [17] "Identity"       "isomap"         "isoMDS"         "MDS"           
    ## [21] "MDS_angle"      "MDS_smacof"     "metaMDS"        "monoMDS"       
    ## [25] "OLO"            "OLO_average"    "OLO_complete"   "OLO_single"    
    ## [29] "OLO_ward"       "optics"         "QAP_2SUM"       "QAP_BAR"       
    ## [33] "QAP_Inertia"    "QAP_LS"         "R2E"            "Random"        
    ## [37] "Reverse"        "Sammon_mapping" "Spectral"       "Spectral_norm" 
    ## [41] "SPIN_NH"        "SPIN_STS"       "TSP"            "VAT"

Details about the method can be found in the [manual page for
`seriate()`](https://mhahsler.r-universe.dev/seriation/doc/manual.html#seriate).

We use a loop to run the function
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md)
with each method and calculate criterion measures which indicate how
good the order is.

``` r

orders <- list()
criterion <- list()
for (m in methods) {
  cat(m)
  tm <- system.time(orders[[m]] <- seriate(d, method = m))
  criterion[[m]] <- data.frame(time = tm[1]+tm[2], rbind(criterion(d, orders[[m]]))) 
  cat(" took", tm[1]+tm[2], "sec.\n")
}
```

    ## ARSA took 0.229 sec.
    ## DendSer took 0.33 sec.
    ## DendSer_ARc took 0.399 sec.
    ## DendSer_BAR took 0.341 sec.
    ## DendSer_LPL took 0.395 sec.
    ## DendSer_PL took 0.395 sec.
    ## GW took 0.213 sec.
    ## GW_average took 0.207 sec.
    ## GW_complete took 0.214 sec.
    ## GW_single took 0.211 sec.
    ## GW_ward took 0.21 sec.
    ## HC took 0.212 sec.
    ## HC_average took 0.212 sec.
    ## HC_complete took 0.212 sec.
    ## HC_single took 0.211 sec.
    ## HC_ward took 0.212 sec.
    ## Identity took 0.211 sec.
    ## isomap

    ##  took 0.245 sec.
    ## isoMDS took 0.211 sec.
    ## MDS took 0.212 sec.
    ## MDS_angle took 0.214 sec.
    ## MDS_smacof took 0.217 sec.
    ## metaMDS took 0.235 sec.
    ## monoMDS took 0.218 sec.
    ## OLO took 0.211 sec.
    ## OLO_average took 0.215 sec.
    ## OLO_complete took 0.217 sec.
    ## OLO_single took 0.216 sec.
    ## OLO_ward took 0.215 sec.
    ## optics took 0.217 sec.
    ## QAP_2SUM took 0.213 sec.
    ## QAP_BAR took 0.217 sec.
    ## QAP_Inertia took 0.218 sec.
    ## QAP_LS took 0.219 sec.
    ## R2E took 0.218 sec.
    ## Random took 0.216 sec.
    ## Reverse took 0.216 sec.
    ## Sammon_mapping took 0.217 sec.
    ## Spectral took 0.218 sec.
    ## Spectral_norm took 0.213 sec.
    ## SPIN_NH took 0.233 sec.
    ## SPIN_STS took 0.229 sec.
    ## TSP took 0.222 sec.
    ## VAT took 0.212 sec.

``` r

criterion <- do.call(rbind, criterion)
```

We align the seriation orders. The reason is that an order 1, 2, 3 and
3, 2, 1 are equivalent and just an artifact of the algorithm. Aligning
will reverse some orders so they are better aligned. Then we sort the
orders from best to worst according to a popular seriation criterion
measure called `Gradient_weighted`.

``` r

orders <- ser_align(orders)
best_to_worse <- order(criterion[["Gradient_weighted"]], decreasing = TRUE)

orders <- orders[best_to_worse]
criterion <- criterion[best_to_worse, ]
```

### Comparison between methods

We can compare the seriation methods by how similar the orders are that
they produce (measured using Spearman). The following code calculates
distances between orders and then performs hierarchical clustering.

``` r

dst <- ser_dist(orders) 
hc <- permute(hclust(dst), order = "OLO", dist = dst)
plot(hc)
```

![](comparison_files/figure-html/unnamed-chunk-6-1.png)

The reordered dendrogram clearly shows a group of methods based on
hierarchical clustering focused on path length and another group that
tries to optimize the other seriation measures.

Here is a table to compare the seriation methods on different criterion
measures. Use the interactive table to sort the methods given different
measures. Note that some are maximized and some should be minimized.
Details about the measures can be found in the [manual page for
`criterion()`](https://mhahsler.r-universe.dev/seriation/doc/manual.html#criterion).

``` r

library(DT)
datatable(round(criterion, 2), extensions = "FixedColumns",
    options = list(paging = TRUE, searching = TRUE, info = FALSE,
      sort = TRUE, scrollX = TRUE, fixedColumns = list(leftColumns = 1))) %>%
    formatRound(columns = colnames(criterion) , mark = "", digits=1)
```

### Visualize the results

Plot the reordered dissimilarity matrices. Dark blocks along the main
diagonal mean that the order reveals a “cluster” of similar objects. The
Iris dataset contains three species, but two of them are very similar,
so we expect to see one smaller block and one larger block.

``` r

for (n in names(orders))
  pimage(d, orders[[n]], main = n , key = FALSE)
```

![](comparison_files/figure-html/unnamed-chunk-8-1.png)![](comparison_files/figure-html/unnamed-chunk-8-2.png)![](comparison_files/figure-html/unnamed-chunk-8-3.png)![](comparison_files/figure-html/unnamed-chunk-8-4.png)![](comparison_files/figure-html/unnamed-chunk-8-5.png)![](comparison_files/figure-html/unnamed-chunk-8-6.png)![](comparison_files/figure-html/unnamed-chunk-8-7.png)![](comparison_files/figure-html/unnamed-chunk-8-8.png)![](comparison_files/figure-html/unnamed-chunk-8-9.png)![](comparison_files/figure-html/unnamed-chunk-8-10.png)![](comparison_files/figure-html/unnamed-chunk-8-11.png)![](comparison_files/figure-html/unnamed-chunk-8-12.png)![](comparison_files/figure-html/unnamed-chunk-8-13.png)![](comparison_files/figure-html/unnamed-chunk-8-14.png)![](comparison_files/figure-html/unnamed-chunk-8-15.png)![](comparison_files/figure-html/unnamed-chunk-8-16.png)![](comparison_files/figure-html/unnamed-chunk-8-17.png)![](comparison_files/figure-html/unnamed-chunk-8-18.png)![](comparison_files/figure-html/unnamed-chunk-8-19.png)![](comparison_files/figure-html/unnamed-chunk-8-20.png)![](comparison_files/figure-html/unnamed-chunk-8-21.png)![](comparison_files/figure-html/unnamed-chunk-8-22.png)![](comparison_files/figure-html/unnamed-chunk-8-23.png)![](comparison_files/figure-html/unnamed-chunk-8-24.png)![](comparison_files/figure-html/unnamed-chunk-8-25.png)![](comparison_files/figure-html/unnamed-chunk-8-26.png)![](comparison_files/figure-html/unnamed-chunk-8-27.png)![](comparison_files/figure-html/unnamed-chunk-8-28.png)![](comparison_files/figure-html/unnamed-chunk-8-29.png)![](comparison_files/figure-html/unnamed-chunk-8-30.png)![](comparison_files/figure-html/unnamed-chunk-8-31.png)![](comparison_files/figure-html/unnamed-chunk-8-32.png)![](comparison_files/figure-html/unnamed-chunk-8-33.png)![](comparison_files/figure-html/unnamed-chunk-8-34.png)![](comparison_files/figure-html/unnamed-chunk-8-35.png)![](comparison_files/figure-html/unnamed-chunk-8-36.png)![](comparison_files/figure-html/unnamed-chunk-8-37.png)![](comparison_files/figure-html/unnamed-chunk-8-38.png)![](comparison_files/figure-html/unnamed-chunk-8-39.png)![](comparison_files/figure-html/unnamed-chunk-8-40.png)![](comparison_files/figure-html/unnamed-chunk-8-41.png)![](comparison_files/figure-html/unnamed-chunk-8-42.png)![](comparison_files/figure-html/unnamed-chunk-8-43.png)![](comparison_files/figure-html/unnamed-chunk-8-44.png)

## Matrix seriation

Matrix seriation reorders rows and columns of a data matrix. We perform
the same steps as for distances in the previous section.

``` r

methods <- sort(list_seriation_methods("matrix"))

# AOE if for correlation matrices only
methods <- setdiff(methods, c("AOE"))
methods 
```

    ##  [1] "BEA"              "BEA_TSP"          "BK_unconstrained" "CA"              
    ##  [5] "Heatmap"          "Identity"         "LLE"              "Mean"            
    ##  [9] "PCA"              "PCA_angle"        "Random"           "Reverse"

Performing seriation.

``` r

orders <- list()
criterion <- list()
for (m in methods) {
  cat(m)
  tm <- system.time(orders[[m]] <- seriate(x, method = m))
  criterion[[m]] <- data.frame(time = tm[1]+tm[2], rbind(criterion(x, orders[[m]])))
  cat(" took", tm[1]+tm[2], "sec.\n")
}
```

    ## BEA took 0.678 sec.
    ## BEA_TSP took 0.661 sec.
    ## BK_unconstrained took 0.226 sec.
    ## CA took 0.224 sec.
    ## Heatmap took 0.661 sec.
    ## Identity took 0.225 sec.
    ## LLE took 0.232 sec.
    ## Mean took 0.224 sec.
    ## PCA took 0.222 sec.
    ## PCA_angle took 0.22 sec.
    ## Random took 0.219 sec.
    ## Reverse took 0.224 sec.

``` r

criterion <- do.call(rbind, criterion)
```

### Comparison between methods

``` r

datatable(round(criterion, 2), extensions = "FixedColumns",
    options = list(paging = TRUE, searching = TRUE, info = FALSE,
      sort = TRUE, scrollX = TRUE, fixedColumns = list(leftColumns = 1))) %>%
    formatRound(columns = colnames(criterion) , mark = "", digits = 1)
```

### Visualize the results

``` r

best_to_worse <- order(criterion[["Moore_stress"]], decreasing = FALSE)

orders <- orders[best_to_worse]
criterion <- criterion[best_to_worse, ]
```

``` r

for (n in names(orders))
  pimage(x, orders[[n]], main = n , key = FALSE)
```

![](comparison_files/figure-html/unnamed-chunk-13-1.png)![](comparison_files/figure-html/unnamed-chunk-13-2.png)![](comparison_files/figure-html/unnamed-chunk-13-3.png)![](comparison_files/figure-html/unnamed-chunk-13-4.png)![](comparison_files/figure-html/unnamed-chunk-13-5.png)![](comparison_files/figure-html/unnamed-chunk-13-6.png)![](comparison_files/figure-html/unnamed-chunk-13-7.png)![](comparison_files/figure-html/unnamed-chunk-13-8.png)![](comparison_files/figure-html/unnamed-chunk-13-9.png)![](comparison_files/figure-html/unnamed-chunk-13-10.png)![](comparison_files/figure-html/unnamed-chunk-13-11.png)![](comparison_files/figure-html/unnamed-chunk-13-12.png)

### 
