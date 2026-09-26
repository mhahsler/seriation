# Seriation Criteria implemented in the R package seriation

This document contains the seriation criteria for judging the quality of
permutations given data implemented in the R package `seriation`.

This list was created for the following version of seriation:

``` r

library("seriation")
packageVersion('seriation')
```

    ## [1] '1.5.8.1'

Register some additional criteria

``` r

register_DendSer()
register_smacof()
```

The criteria are organized by methods that evaluate the permutation
based on distance data or data matrices.

`merit` specified if the criteria function increases with better fit or
if it is formulated as a loss criteria functions (`merit = FALSE`).

If a seriation method directly tries to optimize the criterion, than its
name is specified under “optimized by”. The names can be used as the
`method` argument in
[`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md).

### Criteria for dissimilarity data (dist)

#### 2SUM

2-Sum Criterion: The 2-Sum loss criterion multiplies the similarity
between objects with the squared rank differences (Barnard, Pothen and
Simon, 1993).

- merit: FALSE
- optimized by:
  [QAP_2SUM](http://michael.hahsler.net/seriation/articles/seriation_methods.html#qap_2sum),
  [Spectral](http://michael.hahsler.net/seriation/articles/seriation_methods.html#spectral),
  [Spectral_norm](http://michael.hahsler.net/seriation/articles/seriation_methods.html#spectral_norm)
- additional parameters: no parameters

------------------------------------------------------------------------

#### AR_deviations

Anti-Robinson deviations: The number of violations of the anti-Robinson
form weighted by the deviation (Chen, 2002).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### AR_events

Anti-Robinson events: The number of violations of the anti-Robinson form
(Chen, 2002).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### ARc

Anti-Robinson form cost (Earle and Hurley, 2015).

- merit: FALSE
- optimized by:
  [DendSer_ARc](http://michael.hahsler.net/seriation/articles/seriation_methods.html#dendser_arc)
- registered by: register_DendSer()
- additional parameters: no parameters

------------------------------------------------------------------------

#### BAR

Banded Anti-Robinson form criterion: Measure for closeness to the
anti-Robinson form in a band of size b (Earle and Hurley, 2015).

- merit: FALSE
- optimized by:
  [QAP_BAR](http://michael.hahsler.net/seriation/articles/seriation_methods.html#qap_bar),
  [DendSer_BAR](http://michael.hahsler.net/seriation/articles/seriation_methods.html#dendser_bar)
- additional parameters:

|     | default | help                                     |
|:----|:--------|:-----------------------------------------|
| b   | NULL    | band size defaults to a band of 20% of n |

------------------------------------------------------------------------

#### Gradient_raw

Gradient measure: Evaluates how well distances increase when moving away
from the diagonal of the distance matrix (Hubert et al, 2001).

- merit: TRUE
- optimized by:
  [BBURCG](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bburcg)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Gradient_weighted

Gradient measure (weighted): Evaluates how well distances increase when
moving away from the diagonal of the distance matrix (Hubert et al,
2001).

- merit: TRUE
- optimized by:
  [BBWRCG](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bbwrcg)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Inertia

Inertia criterion: Measures the moment of the inertia of dissimilarity
values around the diagonal of the distance matrix (Caraux and Pinloche,
2005).

- merit: TRUE
- optimized by:
  [QAP_Inertia](http://michael.hahsler.net/seriation/articles/seriation_methods.html#qap_inertia)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Lazy_path_length

Lazy path length: A weighted version of the Hamiltonian path criterion
where later distances are less important (Earl and Hurley, 2015).

- merit: FALSE
- optimized by:
  [DendSer_LPL](http://michael.hahsler.net/seriation/articles/seriation_methods.html#dendser_lpl)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Least_squares

Least squares criterion: The sum of squared differences between
distances and the rank differences (Caraux and Pinloche, 2005).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### LS

Linear Seriation Criterion: Weights the distances with the absolute rank
differences (Hubert and Schultz, 1976).

- merit: TRUE
- optimized by:
  [ARSA](http://michael.hahsler.net/seriation/articles/seriation_methods.html#arsa),
  [QAP_LS](http://michael.hahsler.net/seriation/articles/seriation_methods.html#qap_ls)
- additional parameters: no parameters

------------------------------------------------------------------------

#### MDS_stress

Normalized stress of a configuration given by a seriation order

- merit: FALSE
- optimized by:
  [MDS](http://michael.hahsler.net/seriation/articles/seriation_methods.html#mds),
  [isoMDS](http://michael.hahsler.net/seriation/articles/seriation_methods.html#isomds),
  [Sammon_mapping](http://michael.hahsler.net/seriation/articles/seriation_methods.html#sammon_mapping),
  [monoMDS](http://michael.hahsler.net/seriation/articles/seriation_methods.html#monomds),
  [metaMDS](http://michael.hahsler.net/seriation/articles/seriation_methods.html#metamds)
- additional parameters:

|  | default | help |
|:---|:---|:---|
| warn | FALSE | produce a warning if the 1D MDS fit does not preserve the given order (see ? seriation::uniscale). |

------------------------------------------------------------------------

#### ME

Measure of effectiveness applied to the reordered similarity matrix
(McCormick, 1972).

- merit: TRUE
- optimized by:
  [BEA](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bea),
  [BEA_TSP](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bea_tsp)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Moore_stress

Stress criterion (Moore neighborhood) applied to the reordered
similarity matrix (Niermann, 2005).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### Neumann_stress

Stress criterion (Neumann neighborhood) applied to the reordered
similarity matrix (Niermann, 2005).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### Path_length

Hamiltonian path length: Sum of distances by following the permutation
(Caraux and Pinloche, 2005).

- merit: FALSE
- optimized by:
  [TSP](http://michael.hahsler.net/seriation/articles/seriation_methods.html#tsp),
  [GW](http://michael.hahsler.net/seriation/articles/seriation_methods.html#gw),
  [GW_single](http://michael.hahsler.net/seriation/articles/seriation_methods.html#gw_single),
  [GW_average](http://michael.hahsler.net/seriation/articles/seriation_methods.html#gw_average),
  [GW_complete](http://michael.hahsler.net/seriation/articles/seriation_methods.html#gw_complete),
  [GW_ward](http://michael.hahsler.net/seriation/articles/seriation_methods.html#gw_ward),
  [OLO](http://michael.hahsler.net/seriation/articles/seriation_methods.html#olo),
  [OLO_single](http://michael.hahsler.net/seriation/articles/seriation_methods.html#olo_single),
  [OLO_average](http://michael.hahsler.net/seriation/articles/seriation_methods.html#olo_average),
  [OLO_complete](http://michael.hahsler.net/seriation/articles/seriation_methods.html#olo_complete),
  [OLO_ward](http://michael.hahsler.net/seriation/articles/seriation_methods.html#olo_ward),
  [DendSer_PL](http://michael.hahsler.net/seriation/articles/seriation_methods.html#dendser_pl)
- additional parameters: no parameters

------------------------------------------------------------------------

#### RGAR

Relative generalized anti-Robinson events: Counts Anti-Robinson events
in a variable band of size w around the main diagonal and normalizes by
the maximum of possible events (Tien et al, 2008).

- merit: FALSE
- optimized by: N/A
- additional parameters:

|  | default | help |
|:---|:---|:---|
| w | NULL | window size. Default is to use a pct of 100% of n |
| pct | 100 | specify w as a percentage of n in (0,100\] |
| relative | TRUE | set to FALSE to get the GAR, i.e., the absolute number of AR events in the window. |

------------------------------------------------------------------------

#### Rho

Absolute value of the Spearman rank correlation between original
distances and rank differences of the order.

- merit: TRUE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### smacof_stress0

Stress0 calculated for different transformation types from package
smacof.

- merit: FALSE
- optimized by:
  [MDS_smacof](http://michael.hahsler.net/seriation/articles/seriation_methods.html#mds_smacof)
- registered by: register_smacof()
- additional parameters:

|  | default | help |
|:---|:---|:---|
| type | “ratio” | MDS type (see ? smacof::stress0) |
| warn | FALSE | produce a warning if the 1D MDS fit does not preserve the given order (see ? seriation::uniscale). |
| more | NA | more arguments are passed on to smacof::stress0. |

------------------------------------------------------------------------

### Criteria for data matrices, tables and data.frames

#### Cor_R

Weighted correlation coefficient R: A measure of effectiveness
normalized between -1 and 1 (Deutsch and Martin, 1971).

- merit: TRUE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### ME

Measure of effectiveness (McCormick, 1972).

- merit: TRUE
- optimized by:
  [BEA](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bea),
  [BEA_TSP](http://michael.hahsler.net/seriation/articles/seriation_methods.html#bea_tsp)
- additional parameters: no parameters

------------------------------------------------------------------------

#### Moore_stress

Stress criterion (Moore neighborhood) applied to the reordered matrix
(Niermann, 2005).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

#### Neumann_stress

Stress criterion (Neumann neighborhood) applied to the reordered matrix
(Niermann, 2005).

- merit: FALSE
- optimized by: N/A
- additional parameters: no parameters

------------------------------------------------------------------------

## References

Barnard, S.T., A. Pothen, and H. D. Simon (1993): A Spectral Algorithm
for Envelope Reduction of Sparse Matrices. *In Proceedings of the 1993
ACM/IEEE Conference on Supercomputing,* 493–502. Supercomputing ’93. New
York, NY, USA: ACM.

Caraux, G. and S. Pinloche (2005): Permutmatrix: A Graphical Environment
to Arrange Gene Expression Profiles in Optimal Linear Order,
*Bioinformatics,* **21**(7), 1280–1281.

Chen, C.-H. (2002): Generalized association plots: Information
visualization via iteratively generated correlation matrices,
*Statistica Sinica,* **12**(1), 7–29.

Deutsch, S.B. and J.J. Martin (1971): An ordering algorithm for analysis
of data arrays. *Operational Research,* **19**(6), 1350–1362.
<https://doi.org/10.1287/opre.19.6.1350>

Earle, D. and C.B. Hurley (2015): Advances in Dendrogram Seriation for
Application to Visualization. *Journal of Computational and Graphical
Statistics,* **24**(1), 1–25.
<https://doi.org/10.1080/10618600.2013.874295>

Hahsler, M. (2017): An experimental comparison of seriation methods for
one-mode two-way data. *European Journal of Operational Research,*
**257**, 133–143. <https://doi.org/10.1016/j.ejor.2016.08.066>

Hubert, L. and J. Schultz (1976): Quadratic Assignment as a General Data
Analysis Strategy. *British Journal of Mathematical and Statistical
Psychology,* **29**(2). Blackwell Publishing Ltd. 190–241.
<https://doi.org/10.1111/j.2044-8317.1976.tb00714.x>

Hubert, L., P. Arabie, and J. Meulman (2001): *Combinatorial Data
Analysis: Optimization by Dynamic Programming.* Society for Industrial
Mathematics. <https://doi.org/10.1137/1.9780898718553>

Niermann, S. (2005): Optimizing the Ordering of Tables With Evolutionary
Computation, *The American Statistician,* **59**(1), 41–46.
<https://doi.org/10.1198/000313005X22770>

McCormick, W.T., P.J. Schweitzer and T.W. White (1972): Problem
decomposition and data reorganization by a clustering technique,
*Operations Research,* **20**(5), 993-1009.
<https://doi.org/10.1287/opre.20.5.993>

Robinson, W.S. (1951): A method for chronologically ordering
archaeological deposits, *American Antiquity,* **16**, 293–301.
<https://doi.org/10.2307/276978>

Tien, Y-J., Yun-Shien Lee, Han-Ming Wu and Chun-Houh Chen (2008):
Methods for simultaneously identifying coherent local clusters with
smooth global patterns in gene expression profiles, *BMC
Bioinformatics,* **9**(155), 1–16.
<https://doi.org/10.1186/1471-2105-9-155>
