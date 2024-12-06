# seriation 1.5.7 (12/05/2024)

## New Features
- Added seriation method BK_unconstrained by kbvernon.
- All methods now gracefully handle data with two few objects.
- ser_permutation_vector() now updates method name.

## Bug Fixes
- Fixed label order for seriate.matrix.
- Fixed typo in parameter name for seriation method 
  ARSA (reported by Brian Ripley)

# seriation 1.5.6 (08/19/2024)

## New Features
- Added registered_by field to registries. 

## Changes
- We replaced the FORTRAN implementation for BEA with code from package TSP.
- ME is now calculated using C code.
- optimal.c: updated memory allocation to R allocation.
- stress.c: updated memory allocation to R allocation.

## Bug Fixes
- Added two missing package anchors to palette man page.


# seriation 1.5.5 (04/17/2024)

## Changes
- Updated man pages.

# seriation 1.5.4 (12/11/2023)

## Bug Fixes
- Fixed MDS_angle order for different BLAS implementation giving different
  results for eigen().

# seriation 1.5.3 (11/28/2023)

## New Features
- permute for dendrograms gained parameter dist and accepts now seriation 
  methods.
- Added method "AOE" for correlation matrices.
- registry for seriation methods now contains the name of the seriation criterion
  and a description. seriate_rep now automatically uses the criterion from 
  the registry.
- all seriation methods gained parameter rep.

## Bug Fixes
- optimal.c: use now the correct data type for Rprintf
- Skip deterministic tests on Mac M1 because of numerical differences. 

# seriation 1.5.1 (07/20/2023)

## New Features
- pimage and permute now accept order = TRUE to perform the default seriation.
- hmap gained parameter col_dist to define the color palette used for distance
  matrices.
- hmap dropped parameter showDend and gained parameter plot_margins instead.

## Bug Fixes
- pimage/ggpimage now use zlim correctly to choose the color palette.
- BEA for matrix is now correctly registered as randomized.
- fixed col/row_labels parameter.
- rev() for seriations based on hclust now reverses the dendrogram.
- tests now also accept reverse orders for testing deterministic methods.

# seriation 1.5.0 (07/19/2023)

## New Features
- The seriation registry now contains help information for the seriation method
  parameters.
- New function seriate_best, seriate_rep, and seriate_improve() to easily find a good order for
  randomized algorithms. Parallel execution is supported.
- Seriation method registry has new fields 'randomized' to indicate if an algorithm
  is randomized and can be run several times and 'optimizes' to indicate
  what criterion is optimized. This information is used by seriate_rep.
- Seriation for arrays (including matrix) gained margin parameter.
- tsne and umap can now be used on data matrices.
- get_rank() returns now labels.
- Embedding-based methods now return the order with an attribute called configuration.
- New MDS_stress() function.
- Added register_smacof().
- Added seriation method "Reverse" for dist.
- New seriation methods from vegan: isomap, monoMDS, metaMDS.
- New seriation method "Enumerate" for complete enumeration.
- New seriation method "Mean" for matrix.
- New seriation method "SGD" for distances to improve solutions using stochastic gradient descent.
- New seriation method "LLE" (locally linear embedding) for matrix.
- Heatmap seriation has now special seriation_method "HC_Mean".
- New  seriation criterion "Rho" calculates the absolute Spearman's rank correlation coefficient.
- list_seriation_methods() and list_criterion_methods() gained parameter names_only. 

## Changes
- Seriation methods for MDS are now MDS, isoMDS and Sammon_mapping and have now
  individual control parameters.
- orderplot() is now called plot_config() and can also visualize 2D configurations.
- HC-based seriation: The control parameter method is now linkage so it can be used
  in seriate() in the ...
- Seriation method spectral now also returns the embedding.
- Seriation method simulated annealing is now called "GSA".
- Simplified generics for pimage and ggpimage. Defaults for pimage.dist have changed.
- DendSer methods now return hclust objects.

## Bug Fixes
- fixed labels returned by uniscale()
- FORTRAN: replaced old DFLOAT with DBLE (reported by Brian D. Ripley).

# seriation 1.4.2 (03/07/2023)

## Bug Fixes
-   pimage: ... is now passed on to the seriation method.
-   added missing S3 method registrations.

## New Features
-   methods umap and tsne can now return the embedding.

# seriation 1.4.1 (12/27/2022)

## New Features
-   get_order not consistently returns permutation vectors with names (by david-barnett).

## Bug Fixes
-   criterion.c: replaced enum for bool with <stdbool.h>
-   Additional contributors are not in alphabetical order.

# seriation 1.4.0 (10/21/2022)

## New Features
-   seriate for arrays (including matrices) now returns a complete ser_permutation for all
    dimensions even if margins are specified. For not specified margins, identity permutations
    are returned.
-   added support for tables.
-   added new seriation method CA (correspondence analysis) contributed by Michael Friendly.
-   permute now accepts more than one margin.
-   permute now accepts a seriation method instead of order.

## Bug Fixes
-   seriate.dist now throws correct error upon encountering NAs (by david-barnett)

# seriation 1.3.6 (07/14/2022)

## New Features

-   ggpimage has now a zlim parameter.

## Bug Fixes

-   added register functions back to export (reported by thomasp85).
-   fixed viewports for pimage with colorkey.
-   fixed ggplot diverging color palette direction.


# seriation 1.3.4 (3/16/2022)

## Bug Fixes

-   fixed length calculation in optimal.c

# seriation 1.3.3 (3/3/2022)

## New Features

-   pimage and dissplot gained parameter diag. pimage for dist by
    default does not show the diagonal now.
-   C code now supports long vectors for dist objects.

## User-Visible Changes

-   removed deprecated show functions for the registries.

## Internal Changes

-   we now use roxygen for documentation.
-   added check for long vectors that FORTRAN cannot handle.

# seriation 1.3-2 (2/10/2022)

## Changes

-   improved argument checking for ser_permutation_vector().
-   ggplot uses now standard ggplot2 color palettes.

# seriation 1.3-1 (10/15/2021)

## New Features

-   added seriation based on 1D t-SNE embedding.
-   added seriation based on 1D UMAP embedding.
-   added seriation based on OPTICS.

## User-Visible Changes

-   VAT plots now default to upper_tri = TRUE to show the whole matrix.

# seriation 1.3-0 (06/29/2021)

## User-Visible Changes

-   Plotting
    -   Most plotting functions have now a common interface. This
        changed many parameters.
    -   hmap now uses heatmap from package stats.
    -   dissplot shows now averages in the top triangles.
    -   improved layout (less white space) for grid-based plots.
-   Registry
    -   list_seriation_methods and list_criterion_methods without kind
        return now a list.
    -   show_seriation_methods and show_criterion_methods are deprecated
-   Other Changes
    -   criterion returns now NA with a warning for ME for non-positive
        matrices (used to stop with an error).
    -   dependency dendextend is now only suggested (used for testing).
    -   get_order now returns also labels.
    -   hclust-based seriations now defaults for linkage to complete
        instead of average.

## New Features

-   Plotting
    -   Major refactoring of plotting functions to provide a more
        consistent interface.
    -   added ggplot2-based plots, ggimage, gghmap, ggVAT, ggiVAT,
        ggbertinplot, ggdissplot.
    -   colors are now more consistent and all have bias and power.
-   Seriation methods
    -   seriate for matrix has now method "Heatmap".
    -   seriate now accepts data.frames and used method "heatmap" as the
        default.
    -   added seriation method "Reverse" for reverse identity order.
-   Permutation
    -   permute for matrix-like objects gained parameter margin.
    -   permute for data.frame works now identical to permute for
        matrix.

# seriation 1.2-9 (09/29/2020)

-   removed dependency on methods.
-   added DOIs.

# seriation 1.2-8 (08/27/2019)

## New features

-   get_seriation_method now has better information and also show
    available control parameters.

## Bug Fixes

-   GA: Updated parameter names after change in package ga.

# seriation 1.2-7 (06/07/2019)

## Bug Fixes

-   Added missing void \* to init.c

# seriation 1.2-6 (06/03/2019)

## Bug Fixes

-   Converted print routines in FORTRAN code to dblepr, intpr, etc.
-   seriate.matrix also prints now method name for control verbose =
    TRUE.

# seriation 1.2-5 (05/30/2019)

## Bug Fixes

-   Fixed compilation warnings in FORTRAN code.

# seriation 1.2-4 (05/29/2019)

## New features

-   bertinplot: panel colors can now be specified in highlight and as
    shading.function.

## Bug Fixes

-   bertinplot: fix white squares when frame = TRUE (by Dirk
    Seidensticker).
-   seriation method "BEA" has now a slight code improvement (suggested
    by RichardKav)

# seriation 1.2-3 (02/05/2018)

## Bug Fixes

-   seriation method "BEA" is now not longer masked by "BEA_TSP". Also
    the FORTRAN calls now work.
-   SPIN: making the matrix doubly stochastic now checks all
    rows/columns (reported and fixed by cerebis)

# seriation 1.2-2 (05/08/2017)

## New features

-   Added new seriation method SA which provides simulated annealing for
    all criterion measures.
-   Added criterion Cor_R (ME for the moment ordering algorithm by
    Deutsch and Martin).
-   Added uniscale to produce a unidimensional scaling configuration
    given a distance matrix and a permutation.
-   Criterion gained parameter force_loss (default is FALSE). Merit
    measures are converted into loss values by multiplying with -1.
-   Added Supreme Court dataset.

## Changes and Bug Fixes

-   Default for seriate (dist) and dissplot is now "Spectral" since it
    gives a better tradeoff between quality and speed.
-   Seriation method ARSA's control argument nreps is now for
    consistency called reps.
-   Criterion: dist objects are now automatically converted into a
    similarity matrix for ME, Moore_stress and Neumann_stress.
-   pimage now suppresses the color key for logical matrices and checks
    for all NAs and infinite entries.
-   Correction: ARSA minimizes the linear seriation criterion (man page
    and vignette).

# seriation 1.2-1 (08/06/2016)

## New features

-   Added new distance measure called absolute pairwise rank
    differences.

## Changes and Bug Fixes

-   The default setting for ser_dist and ser_cor is now reverse is TRUE.
-   pimage does now work with matrices containing only a single value.
-   control parameters for method TSP are now correctly passed on
    (reported by David Aliyev).

# seriation 1.2-0 (2/22/2016)

## New features

-   RGAR gained parameter pct to specify the window as a percentage.
-   Added the lazy path length criterion.
-   Added the banded anti-Robinson form (BAR) criterion.
-   Added QAP_Inertia and QAP_BAR solver.
-   Added DendSer using register_DendSer().
-   Added GA using register_GA().

## Changes and Bug Fixes

-   Fixed RGAR (w needs to be in [2,n-1]).
-   Registry now warns and modifies entries with the same name.
-   Registry now lists methods in alphabetical order.
-   Seriation method alias Chen was removed. Use R2E.

# seriation 1.1-3 (12/18/2015)

-   Added is.robinson to recognize (pre) Robinson matrices.
-   Added random.robinson to create random Robinson matrices.
-   Added seriation methods "QAP_LS" and "QAP_2SUM" (QAP-based
    seriation).
-   Added criteria "LS" and "2SUM" from QAP-based seriation.
-   Fixed Spectral_norm seriation.
-   hmap now honors zlim also in dendrogram-based maps.
-   hmap gained option sym for seriation based maps. showdist can now be
    one of "none" (default), "row", "column", or "both".
-   ser_cor and ser_dist gained parameter y. ser_cor gained parameter
    test to perform tests for association.
-   Added permute method for hclust and dendrogram objects.

# seriation 1.1-2 (8/23/2015)

-   Argument (control and ...) check warns now instead of throwing an
    error.
-   seriation_dist, seriation_cor and seriation_align are now shortened
    to ser_dist, ser_cor and ser_align.
-   Method "ppc" is now faster and also available in ser_cor.
-   Fixed ser_cor for "spearman" and "Kendall" (uses now rank
    correctly).
-   ser_cor and ser_dist gained parameter reverse to indicate that
    permuations are also tried in reverse and the best value is
    reported.

# seriation 1.1-1 (7/1/2015)

-   get_permutation_matrix added.
-   seriation_dist measure "ppc" (positional proximity coefficient)
    added.
-   Fixed bug with permute and ser_permutation_vectors.
-   Identity permutations (NA) give now an error for get_order and
    get_permutation_matrix.
-   Fixed imports for non-base R packages.

# seriation 1.1-0 (06/09/2015)

-   Seriation method 'Identity' added.

-   Seriation method 'Random' added.

-   Seriation method 'VAT' added.

-   Seriation methods 'Spectral' and 'Spectral_norm' added.

-   Seriation methods 'PCA_angle' and 'MDS_angle' added.

-   Seriation methods 'SPIN_NH' and 'SPIN_STS' added.

-   Several aliases for seriation methods added.

-   Criterion 'RWGAR' added.

-   permutation_matrix2vector and permutation_vector2matrix added.

-   Identity permutation (value NA) added.

-   ser_permutation and ser_permutation_vector can now be used
    interchangeably,

-   get_rank for permutation vectors added.

-   seriation_dist and seriation_alignment to calculate dissimilarities
    between seriation orders added.

-   Wood data set added.

-   

    # Chameleon data sets added.

-   create_lines_data, create_ordered_data added.

-   pimage, hmap and dissplot: Simplified and made interfaces more
    consistent (all use now zlim, consistent default color palettes).

-   pimage gained axes and prop; NA in matrix now works.

-   seriation checks now control arguments consistently.

-   We use now package registry to manage methods.

-   reorder for hclust added.

-   iVAT with path distance added.

-   color palettes (bluered, greenred, grays) added.

-   Improved speed of C code.

-   Problem with testthat file names fixed.

-   bburg.f/bbwrg.f: memory access problem fixed.

# seriation 1.0-14 (12/02/2014)

-   arsa.f: removed 0 flag in rand() so it compiles under AIX (reported
    by Lei Zhang)
-   arsa.f/bburg.f/bbwrg.f: calls now R RNG to be compatible with
    certain compilers (e.g., Intel FORTRAN) (reported by Rohan Shah)

# seriation 1.0-13 (3/11/2014)

-   Fixed dependence on MASS

# seriation 1.0-12 (2/18/2014)

-   ser_permutation_vectors can now be reversed with rev
-   get_order: removed the weird labels.
-   we use now testthat
-   fixed bug with intra-cluster ordering using silhouette width
    (reported by Bettina Gruen)
-   Cleaned up dependencies: TSP, grid, cluster, gclus and colorspace
    are now imports instead of dependencies.

# seriation 1.0-11 (9/6/2013)

-   service release.

# seriation 1.0-10 (2/15/2013)

-   pimage has now a colorkey and a range argument
-   fixed bug in ARSA when the distance matrix contains all 0s
-   added PACKAGE argument to .Fortran calls

# seriation 1.0-8 and 1.0-9 (11/6/2012)

-   get_order: labels are now in the correct order (Bug report by Crt
    Ahlin)
-   Replaced FORTRAN I/O with R I/O for verb=TRUE
-   Fixed pop/newpage bug in pimage.dist (reported by Bettina Gruen)

# seriation 1.0-7 (9/25/2012)

-   Fixed out-of-bounds bug in arsa.f (reported by Rohan Shah)
-   Fixed out-of-bounds bug in bburcg.f

# seriation 1.0-6 (10/19/2011)

-   removed deprecated parameter gamma for dissplot()

# seriation 1.0-5 (9/2/2011)

-   bertinplot(): fixed representation for 0, neg. values and highlight.
    (Bug report by G. Sawitzki).
-   bertinplot(): added panel.blocks and option for shading
-   bertinplot(): added bertin_cut_line()

# seriation 1.0-4 (6/28/2011)

-   pimage() now uses grid.raster.
-   dissplot() now uses grid.raster.

# seriation 1.0-3 (1/14/2011)

-   improved validity check for permutations and added check for dist
    with neg. entries to seriate.dist.

# seriation 1.0-2 (3/13/2010)

-   service release

# seriation 1.0-1 (8/25/2009)

-   added drop=FALSE in permute for matrix.
-   fixed reordering for labels.
-   added permute for character.
-   added different methods to calculate between cluster dissimilarities
    (min, max, avg, Hausdorff).
-   dissplot has now additional options hue, power, gamma, flip and
    changed behavior for averages. dissplot depends now on colorspace.

# Version 1.0-0 (3/24/2009)

-   many changes and first stable release.

# Version 0.1-1 (9/1/2007)

-   Initial beta release.
