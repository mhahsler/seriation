# Package index

## Seriation

Create a seriation order which is a permutation of objects.

- [`seriate()`](http://michael.hahsler.net/seriation/reference/seriate.md)
  : Seriate Dissimilarity Matrices, Matrices or Arrays
- [`seriate_best()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)
  [`seriate_rep()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)
  [`seriate_improve()`](http://michael.hahsler.net/seriation/reference/seriate_best.md)
  : Best Seriation
- [`register_DendSer()`](http://michael.hahsler.net/seriation/reference/register_DendSer.md)
  : Register Seriation Methods from Package DendSer
- [`register_GA()`](http://michael.hahsler.net/seriation/reference/register_GA.md)
  [`gaperm_mixedMutation()`](http://michael.hahsler.net/seriation/reference/register_GA.md)
  : Register a Genetic Algorithm Seriation Method
- [`register_optics()`](http://michael.hahsler.net/seriation/reference/register_optics.md)
  : Register Seriation Based on OPTICS
- [`register_smacof()`](http://michael.hahsler.net/seriation/reference/register_smacof.md)
  : Register Seriation Methods from Package smacof
- [`register_tsne()`](http://michael.hahsler.net/seriation/reference/register_tsne.md)
  : Register Seriation Based on 1D t-SNE
- [`register_umap()`](http://michael.hahsler.net/seriation/reference/register_umap.md)
  : Register Seriation Based on 1D UMAP
- [`registry_seriate`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`list_seriation_methods()`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`get_seriation_method()`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`set_seriation_method()`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`print(`*`<seriation_method>`*`)`](http://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  : Registry for Seriation Methods

## Seriation Criterion

Calculate seriation creteria for a seriation order.

- [`criterion()`](http://michael.hahsler.net/seriation/reference/criterion.md)
  : Criterion for a Loss/Merit Function for Data Given a Permutation
- [`registry_criterion`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`list_criterion_methods()`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`get_criterion_method()`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`set_criterion_method()`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`print(`*`<criterion_method>`*`)`](http://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  : Registry for Criterion Methods

## Permutations

Work with permutations returned by seriation methods.

- [`get_order()`](http://michael.hahsler.net/seriation/reference/get_order.md)
  [`get_rank()`](http://michael.hahsler.net/seriation/reference/get_order.md)
  [`get_permutation_matrix()`](http://michael.hahsler.net/seriation/reference/get_order.md)
  : Extracting Order Information from a Permutation Object
- [`permute()`](http://michael.hahsler.net/seriation/reference/permute.md)
  : Permute the Order in Various Objects
- [`permutation_vector2matrix()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md)
  [`permutation_matrix2vector()`](http://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md)
  : Conversion Between Permutation Vector and Permutation Matrix
- [`reorder(`*`<hclust>`*`)`](http://michael.hahsler.net/seriation/reference/reorder.hclust.md)
  : Reorder Dendrograms using Optimal Leaf Ordering
- [`ser_dist()`](http://michael.hahsler.net/seriation/reference/ser_dist.md)
  [`ser_cor()`](http://michael.hahsler.net/seriation/reference/ser_dist.md)
  [`ser_align()`](http://michael.hahsler.net/seriation/reference/ser_dist.md)
  : Dissimilarities and Correlations Between Seriation Orders
- [`ser_permutation()`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`print(`*`<ser_permutation>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`summary(`*`<ser_permutation>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`c(`*`<ser_permutation>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`` `[`( ``*`<ser_permutation>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation.md)
  : Class ser_permutation – A Collection of Permutation Vectors for
  Seriation
- [`ser_permutation_vector()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`c(`*`<ser_permutation_vector>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`rev(`*`<ser_permutation_vector>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`get_method()`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`length(`*`<ser_permutation_vector>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`print(`*`<ser_permutation_vector>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`summary(`*`<ser_permutation_vector>`*`)`](http://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  : Class ser_permutation_vector – A Single Permutation Vector for
  Seriation

## Visualizations

Visualize seriated data to discover patterns.

- [`VAT()`](http://michael.hahsler.net/seriation/reference/VAT.md)
  [`iVAT()`](http://michael.hahsler.net/seriation/reference/VAT.md)
  [`path_dist()`](http://michael.hahsler.net/seriation/reference/VAT.md)
  [`ggVAT()`](http://michael.hahsler.net/seriation/reference/VAT.md)
  [`ggiVAT()`](http://michael.hahsler.net/seriation/reference/VAT.md) :
  Visual Analysis for Cluster Tendency Assessment (VAT/iVAT)
- [`bertinplot()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.bars()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.circles()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.rectangles()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.squares()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.tiles()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.blocks()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.lines()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`bertin_cut_line()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`ggbertinplot()`](http://michael.hahsler.net/seriation/reference/bertinplot.md)
  : Plot a Bertin Matrix
- [`dissplot()`](http://michael.hahsler.net/seriation/reference/dissplot.md)
  [`plot(`*`<reordered_cluster_dissimilarity_matrix>`*`)`](http://michael.hahsler.net/seriation/reference/dissplot.md)
  [`print(`*`<reordered_cluster_dissimilarity_matrix>`*`)`](http://michael.hahsler.net/seriation/reference/dissplot.md)
  [`ggdissplot()`](http://michael.hahsler.net/seriation/reference/dissplot.md)
  : Dissimilarity Plot
- [`hmap()`](http://michael.hahsler.net/seriation/reference/hmap.md)
  [`gghmap()`](http://michael.hahsler.net/seriation/reference/hmap.md) :
  Plot Heat Map Reordered Using Seriation
- [`bluered()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`greenred()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`reds()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`blues()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`greens()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`greys()`](http://michael.hahsler.net/seriation/reference/palette.md)
  [`grays()`](http://michael.hahsler.net/seriation/reference/palette.md)
  : Different Useful Color Palettes
- [`pimage()`](http://michael.hahsler.net/seriation/reference/pimage.md)
  [`ggpimage()`](http://michael.hahsler.net/seriation/reference/pimage.md)
  : Permutation Image Plot

## Data Sets

- [`Chameleon`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds4`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds5`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds7`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds8`](http://michael.hahsler.net/seriation/reference/Chameleon.md)
  : 2D Data Sets used for the CHAMELEON Clustering Algorithm
- [`Irish`](http://michael.hahsler.net/seriation/reference/Irish.md) :
  Irish Referendum Data Set
- [`Munsingen`](http://michael.hahsler.net/seriation/reference/Munsingen.md)
  : Hodson's Munsingen Data Set
- [`Psych24`](http://michael.hahsler.net/seriation/reference/Psych24.md)
  : Results of 24 Psychological Test for 8th Grade Students
- [`SupremeCourt`](http://michael.hahsler.net/seriation/reference/SupremeCourt.md)
  : Voting Patterns in the Second Rehnquist U.S. Supreme Court
- [`Townships`](http://michael.hahsler.net/seriation/reference/Townships.md)
  : Bertin's Characteristics of Townships
- [`Wood`](http://michael.hahsler.net/seriation/reference/Wood.md) :
  Gene Expression Data for Wood Formation in Poplar Trees
- [`Zoo`](http://michael.hahsler.net/seriation/reference/Zoo.md) : Zoo
  Data Set
- [`create_lines_data()`](http://michael.hahsler.net/seriation/reference/create_lines_data.md)
  [`create_ordered_data()`](http://michael.hahsler.net/seriation/reference/create_lines_data.md)
  : Create Simulated Data for Seriation Evaluation
- [`is.robinson()`](http://michael.hahsler.net/seriation/reference/is.robinson.md)
  [`random.robinson()`](http://michael.hahsler.net/seriation/reference/is.robinson.md)
  : Create and Recognize Robinson and Pre-Robinson Matrices

## Helper

- [`LS_swap()`](http://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_insert()`](http://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_reverse()`](http://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_mixed()`](http://michael.hahsler.net/seriation/reference/LS.md) :
  Neighborhood functions for Seriation Method SA
- [`lle()`](http://michael.hahsler.net/seriation/reference/lle.md) :
  Locally Linear Embedding (LLE)
- [`uniscale()`](http://michael.hahsler.net/seriation/reference/uniscale.md)
  [`MDS_stress()`](http://michael.hahsler.net/seriation/reference/uniscale.md)
  [`get_config()`](http://michael.hahsler.net/seriation/reference/uniscale.md)
  [`plot_config()`](http://michael.hahsler.net/seriation/reference/uniscale.md)
  : Fit an Unidimensional Scaling for a Seriation Order
