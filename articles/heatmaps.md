# Heatmaps with Package Seriation

## Introduction

A [Heatmap](https://en.wikipedia.org/wiki/Heat_map) uses colored tiles
to represent the values in a data matrix. Patterns can be easily seen if
the rows and columns are appropriately reordered. There are many ways to
reorder a matrix, and the order has a significant impact on the
visualization’s usefulness. The package `seriation` implements a large
number of reordering methods (see: the [list with all implemented
seriation
methods](https://mhahsler.github.io/seriation/seriation_methods.html)).
`seriation` also provides a set of functions to display reordered
heatmaps:

- [`pimage()`](http://michael.hahsler.net/seriation/reference/pimage.md)
- [`hmap()`](http://michael.hahsler.net/seriation/reference/hmap.md)
- [`gghmap()`](http://michael.hahsler.net/seriation/reference/hmap.md)

How to cite the `seriation` package:

> Hahsler M, Hornik K, Buchta C (2008). “Getting things in order: An
> introduction to the R package seriation.” *Journal of Statistical
> Software*, *25*(3), 1-34. ISSN 1548-7660, <doi:10.18637/jss.v025.i03>
> <https://doi.org/10.18637/jss.v025.i03>.

### Prepare the data

As an example, we use the `Wood` dataset with the normalized gene
expression data (a sample of 136 genes) for wood formation in poplar
trees in 6 locations. In case the data already has some order, we
randomly reorder rows and columns for this example.

``` r

if (!require("seriation")) install.packages("seriation")
```

    ## Loading required package: seriation

``` r

library("seriation")
data("Wood")
Wood <- Wood[sample(nrow(Wood)), sample(ncol(Wood))]
dim(Wood)
```

    ## [1] 136   6

``` r

DT::datatable(round(Wood, 2))
```

Here is a simple heatmap without reordering. No structure is visible.

``` r

pimage(Wood)
```

![](heatmaps_files/figure-html/unnamed-chunk-2-1.png)

### Reordering in seriation

Many seriation methods are available. The [manual page for
seriate()](https://mhahsler.r-universe.dev/seriation/doc/manual.html#seriate)
describes the methods available in package `seriation`.

Methods of interest for heatmaps are dendrogram leaf order-based methods
applied to rows and columns. This is done using `method = "heatmap"`.
The actual seriation method can be passed on as parameter
`seriation_method`, but it has a suitable default if it is omitted. Here
is an example:

``` r

o <- seriate(Wood, method = "Heatmap", seriation_method = "HC_Mean")
o
```

    ## object of class 'ser_permutation', 'list'
    ## contains permutation vectors for 2-mode data
    ## 
    ##   vector length seriation method
    ## 1           136          Heatmap
    ## 2             6          Heatmap

This is the order for rows and columns. The method `heatmap`
automatically performs hierarchical clustering and then applies the
seriation method for reordering of dendrogram leaves. Here we use the
row/column mean to reorder the dendrogram. The resulting order (2 means
second dimension, i.e., columns) can be shown, and the reordered
dendrogram and a reordered image can be plotted.

``` r

get_order(o, 2)
```

    ## B A P E C D 
    ## 2 3 5 4 1 6

``` r

plot(o[[2]])
pimage(Wood, order = o)
```

![](heatmaps_files/figure-html/unnamed-chunk-5-1.png)![](heatmaps_files/figure-html/unnamed-chunk-5-2.png)

## Built-in heatmap function

Package `seriation` has several functions to display heatmaps.

### Without dendrograms: pimage

The permutation image plot in `seriation` provides a simple heatmap. The
order argument not only accepts a seriation order, but also a seriation
method.

``` r

pimage(Wood, order = "Heatmap", seriation_method = "HC_complete", 
       main = "Wood (hierarchical clustering)")
pimage(Wood, order = "Heatmap", seriation_method = "HC_Mean", 
       main = "Wood (reorder by row/col mean)")
pimage(Wood, order = "Heatmap", seriation_method = "GW_complete", 
       main = "Wood (reorder by Gruvaeus and Wainer heuristic)")
```

    ## Registered S3 method overwritten by 'gclus':
    ##   method         from     
    ##   reorder.hclust seriation

``` r

pimage(Wood, order = "Heatmap", 
       main = "Wood (default - optimal leaf ordering)")
```

![](heatmaps_files/figure-html/unnamed-chunk-6-1.png)![](heatmaps_files/figure-html/unnamed-chunk-6-2.png)![](heatmaps_files/figure-html/unnamed-chunk-6-3.png)![](heatmaps_files/figure-html/unnamed-chunk-6-4.png)

### With dendrograms: hmap

Here are some typical reordering schemes.

``` r

hmap(Wood, method = "HC_complete", main = "Wood (hierarchical clustering)")
hmap(Wood, method = "HC_Mean", main = "Wood (reorder by row/col mean)")
hmap(Wood, method = "GW_complete", main = "Wood (reorder by Gruvaeus and Wainer heuristic)")
hmap(Wood, method = "OLO_complete", main = "Wood (opt. leaf ordering)")
```

![](heatmaps_files/figure-html/unnamed-chunk-7-1.png)![](heatmaps_files/figure-html/unnamed-chunk-7-2.png)![](heatmaps_files/figure-html/unnamed-chunk-7-3.png)![](heatmaps_files/figure-html/unnamed-chunk-7-4.png)

Different linkage types can be added in the method name.

Package `DendSer` offers more dendrogram seriation methods. These
methods can be registered using \`register_DendSer()\`\`

``` r

register_DendSer()
```

    ## Registering new seriation method 'DendSer' for 'dist'register_DendSer()

    ## Registering new seriation method 'DendSer_BAR' for 'dist'register_DendSer()

    ## Registering new seriation method 'DendSer_PL' for 'dist'register_DendSer()

    ## Registering new seriation method 'DendSer_LPL' for 'dist'register_DendSer()

    ## Registering new seriation method 'DendSer_ARc' for 'dist'register_DendSer()

    ## Registering new seriation criterion 'ARc' for 'dist' using register_DendSer()

``` r

hmap(Wood, method = "DendSer_BAR", main = "Wood (banded anti-Robinson)")
```

![](heatmaps_files/figure-html/unnamed-chunk-8-1.png)

### With distance matrices instead of dendrograms: hmap

Instead of dendrograms, also reordered distance matrices can be
displayed. Dark block around the diagonal indicate the cluster
structure.

``` r

hmap(Wood, method = "HC_complete", 
     plot_margins = "distances",
     main = "Wood (hierarchical clustering)")
```

![](heatmaps_files/figure-html/unnamed-chunk-9-1.png)

Also non-dendrogram-based reordering methods can be used. These methods
reorder rows and columns. Instead of the dendrograms, reordered distance
matrices are shown.

``` r

hmap(Wood, method = "MDS", main = "Wood (MDS)")
hmap(Wood, method = "MDS_angle", main = "Wood (Angle in 2D MDS space)")
hmap(Wood, method = "R2E", main = "Wood (Rank 2 ellipse seriation)")
hmap(Wood, method = "TSP", main = "Wood (Traveling salesperson)")
```

![](heatmaps_files/figure-html/unnamed-chunk-10-1.png)![](heatmaps_files/figure-html/unnamed-chunk-10-2.png)![](heatmaps_files/figure-html/unnamed-chunk-10-3.png)![](heatmaps_files/figure-html/unnamed-chunk-10-4.png)

### colors with pimage and hmap

``` r

hmap(Wood, col = grays())
hmap(Wood, col = greenred())

hmap(Wood, col = colorRampPalette(c("brown", "orange", "red"))( 100 ) )

if (!require("viridis")) install.packages("viridis")
```

    ## Loading required package: viridis

    ## Loading required package: viridisLite

``` r

hmap(Wood, col = viridis::viridis(100))
```

![](heatmaps_files/figure-html/unnamed-chunk-11-1.png)![](heatmaps_files/figure-html/unnamed-chunk-11-2.png)![](heatmaps_files/figure-html/unnamed-chunk-11-3.png)![](heatmaps_files/figure-html/unnamed-chunk-11-4.png)

There are many other packages to create color palettes in R like
`RColorBrewer` or `colorspaces`.

### ggplot2

All options are also available for `ggplot2` using
[`gghmap()`](http://michael.hahsler.net/seriation/reference/hmap.md).
Currently there is no support to display dendrograms.

``` r

if (!require("ggplot2")) install.packages("ggplot2")
```

    ## Loading required package: ggplot2

``` r

library(ggplot2)
gghmap(Wood, method = "OLO")
```

![](heatmaps_files/figure-html/unnamed-chunk-12-1.png)

## Using seriation with other packages

The package `seriation` can be used to compute reordering for other
heatmap packages.

### heatmap in package stats

This is R’s standard heatmap function.
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md)
can be used to supply the reordered dendrogram.

``` r

o <- seriate(Wood, method = "Heatmap", seriation_method = "OLO")
heatmap(Wood, Rowv = as.dendrogram(o[[1]]), Colv = as.dendrogram(o[[2]]))
```

![](heatmaps_files/figure-html/unnamed-chunk-13-1.png)

We can also supply the rank order for any seriation method as weights
and the dendrogram will be reordered as close as possible to the
seriation order.

``` r

o <- seriate(Wood, method = "Heatmap", seriation_method = "Spectral")
heatmap(Wood, Rowv =  get_rank(o, 1), Colv =  get_rank(o, 2))
```

![](heatmaps_files/figure-html/unnamed-chunk-14-1.png)

### Package heatmaply

The package creates interactive heatmaps. It already uses the package
`seriation` for parameter `seriate` and supports the `OLO` and `GW`
methods.

``` r

if (!suppressMessages(require("heatmaply"))) install.packages("heatmaply")
```

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'heatmaply'

    ## Installing package into '/home/runner/work/_temp/Library'
    ## (as 'lib' is unspecified)

    ## also installing the dependencies 'httr', 'plyr', 'plotly', 'reshape2', 'webshot', 'assertthat', 'egg'

``` r

library("heatmaply")
```

    ## Loading required package: plotly

    ## 
    ## Attaching package: 'plotly'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     last_plot

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following object is masked from 'package:graphics':
    ## 
    ##     layout

    ## 
    ## ======================
    ## Welcome to heatmaply version 1.6.0
    ## 
    ## Type citation('heatmaply') for how to cite the package.
    ## Type ?heatmaply for the main documentation.
    ## 
    ## The github page is: https://github.com/talgalili/heatmaply/
    ## Please submit your suggestions and bug-reports at: https://github.com/talgalili/heatmaply/issues
    ## You may ask questions at stackoverflow, use the r and heatmaply tags: 
    ##   https://stackoverflow.com/questions/tagged/heatmaply
    ## ======================

``` r

heatmaply(Wood, seriate = "none", main = "HC")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the dendextend package.
    ##   Please report the issue at <https://github.com/talgalili/dendextend/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r

heatmaply(Wood, seriate = "OLO", main = "OLO")
```

Any dendrogram-based seriation method from `seriation` can be supplied.

``` r

o <- seriate(Wood, method = "Heatmap", seriation_method = "OLO_ward")
heatmaply(Wood, Rowv = o[[1]], Colv = o[[2]], main = "OLO (Ward)")
```

``` r

o <- seriate(Wood, method = "Heatmap", seriation_method = "Spectral")
heatmaply(Wood, Rowv = get_rank(o, 1), Colv = get_rank(o, 2), main = "Spectral")
```
