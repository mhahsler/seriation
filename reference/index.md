# Package index

## Seriation Methods

Create a seriation order which is a permutation of objects.

- [`seriate()`](https://michael.hahsler.net/seriation/reference/seriate.md)
  : Seriate Dissimilarity Matrices, Matrices or Arrays
- [`seriate_best()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)
  [`seriate_rep()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)
  [`seriate_improve()`](https://michael.hahsler.net/seriation/reference/seriate_best.md)
  : Best Seriation
- [`register_DendSer()`](https://michael.hahsler.net/seriation/reference/register_DendSer.md)
  : Register Seriation Methods from Package DendSer
- [`register_GA()`](https://michael.hahsler.net/seriation/reference/register_GA.md)
  [`gaperm_mixedMutation()`](https://michael.hahsler.net/seriation/reference/register_GA.md)
  : Register a Genetic Algorithm Seriation Method
- [`register_optics()`](https://michael.hahsler.net/seriation/reference/register_optics.md)
  : Register Seriation Based on OPTICS
- [`register_smacof()`](https://michael.hahsler.net/seriation/reference/register_smacof.md)
  : Register Seriation Methods from Package smacof
- [`register_tsne()`](https://michael.hahsler.net/seriation/reference/register_tsne.md)
  : Register Seriation Based on 1D t-SNE
- [`register_umap()`](https://michael.hahsler.net/seriation/reference/register_umap.md)
  : Register Seriation Based on 1D UMAP
- [`register_vegan()`](https://michael.hahsler.net/seriation/reference/register_vegan.md)
  : Register Seriation Methods from Package vegan
- [`registry_seriate`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`list_seriation_methods()`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`get_seriation_method()`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`set_seriation_method()`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  [`print(`*`<seriation_method>`*`)`](https://michael.hahsler.net/seriation/reference/registry_for_seriation_methods.md)
  : Registry for Seriation Methods

## Seriation Criterion

Calculate seriation creteria for a seriation order.

- [`criterion()`](https://michael.hahsler.net/seriation/reference/criterion.md)
  : Criterion for a Loss/Merit Function for Data Given a Permutation
- [`registry_criterion`](https://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`list_criterion_methods()`](https://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`get_criterion_method()`](https://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`set_criterion_method()`](https://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  [`print(`*`<criterion_method>`*`)`](https://michael.hahsler.net/seriation/reference/registry_for_criterion_methods.md)
  : Registry for Criterion Methods

## Permutations

Work with permutations returned by seriation methods.

- [`get_order()`](https://michael.hahsler.net/seriation/reference/get_order.md)
  [`get_rank()`](https://michael.hahsler.net/seriation/reference/get_order.md)
  [`get_permutation_matrix()`](https://michael.hahsler.net/seriation/reference/get_order.md)
  : Extracting Order Information from a Permutation Object
- [`permute()`](https://michael.hahsler.net/seriation/reference/permute.md)
  : Permute the Order in Various Objects
- [`permutation_vector2matrix()`](https://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md)
  [`permutation_matrix2vector()`](https://michael.hahsler.net/seriation/reference/permutation_vector2matrix.md)
  : Conversion Between Permutation Vector and Permutation Matrix
- [`reorder(`*`<hclust>`*`)`](https://michael.hahsler.net/seriation/reference/reorder.hclust.md)
  : Reorder Dendrograms using Optimal Leaf Ordering
- [`ser_dist()`](https://michael.hahsler.net/seriation/reference/ser_dist.md)
  [`ser_cor()`](https://michael.hahsler.net/seriation/reference/ser_dist.md)
  [`ser_align()`](https://michael.hahsler.net/seriation/reference/ser_dist.md)
  : Dissimilarities and Correlations Between Seriation Orders
- [`ser_permutation()`](https://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`print(`*`<ser_permutation>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`summary(`*`<ser_permutation>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`c(`*`<ser_permutation>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation.md)
  [`` `[`( ``*`<ser_permutation>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation.md)
  : Class ser_permutation – A Collection of Permutation Vectors for
  Seriation
- [`ser_permutation_vector()`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`c(`*`<ser_permutation_vector>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`rev(`*`<ser_permutation_vector>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`get_method()`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`length(`*`<ser_permutation_vector>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`print(`*`<ser_permutation_vector>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  [`summary(`*`<ser_permutation_vector>`*`)`](https://michael.hahsler.net/seriation/reference/ser_permutation_vector.md)
  : Class ser_permutation_vector – A Single Permutation Vector for
  Seriation

## Visualizations

Visualize seriated data to discover patterns.

- [`VAT()`](https://michael.hahsler.net/seriation/reference/VAT.md)
  [`iVAT()`](https://michael.hahsler.net/seriation/reference/VAT.md)
  [`path_dist()`](https://michael.hahsler.net/seriation/reference/VAT.md)
  [`ggVAT()`](https://michael.hahsler.net/seriation/reference/VAT.md)
  [`ggiVAT()`](https://michael.hahsler.net/seriation/reference/VAT.md) :
  Visual Analysis for Cluster Tendency Assessment (VAT/iVAT)
- [`bertinplot()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.bars()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.circles()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.rectangles()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.squares()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.tiles()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.blocks()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`panel.lines()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`bertin_cut_line()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  [`ggbertinplot()`](https://michael.hahsler.net/seriation/reference/bertinplot.md)
  : Plot a Bertin Matrix
- [`dissplot()`](https://michael.hahsler.net/seriation/reference/dissplot.md)
  [`plot(`*`<reordered_cluster_dissimilarity_matrix>`*`)`](https://michael.hahsler.net/seriation/reference/dissplot.md)
  [`print(`*`<reordered_cluster_dissimilarity_matrix>`*`)`](https://michael.hahsler.net/seriation/reference/dissplot.md)
  [`ggdissplot()`](https://michael.hahsler.net/seriation/reference/dissplot.md)
  : Dissimilarity Plot
- [`hmap()`](https://michael.hahsler.net/seriation/reference/hmap.md)
  [`gghmap()`](https://michael.hahsler.net/seriation/reference/hmap.md)
  : Plot Heat Map Reordered Using Seriation
- [`bluered()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`greenred()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`reds()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`blues()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`greens()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`greys()`](https://michael.hahsler.net/seriation/reference/palette.md)
  [`grays()`](https://michael.hahsler.net/seriation/reference/palette.md)
  : Different Useful Color Palettes
- [`pimage()`](https://michael.hahsler.net/seriation/reference/pimage.md)
  [`ggpimage()`](https://michael.hahsler.net/seriation/reference/pimage.md)
  : Permutation Image Plot

## Data Sets

- [`Chameleon`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds4`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds5`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds7`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  [`chameleon_ds8`](https://michael.hahsler.net/seriation/reference/Chameleon.md)
  : 2D Data Sets used for the CHAMELEON Clustering Algorithm
- [`Irish`](https://michael.hahsler.net/seriation/reference/Irish.md) :
  Irish Referendum Data Set
- [`Munsingen`](https://michael.hahsler.net/seriation/reference/Munsingen.md)
  : Hodson's Munsingen Data Set
- [`Psych24`](https://michael.hahsler.net/seriation/reference/Psych24.md)
  : Results of 24 Psychological Test for 8th Grade Students
- [`SupremeCourt`](https://michael.hahsler.net/seriation/reference/SupremeCourt.md)
  : Voting Patterns in the Second Rehnquist U.S. Supreme Court
- [`Townships`](https://michael.hahsler.net/seriation/reference/Townships.md)
  : Bertin's Characteristics of Townships
- [`Wood`](https://michael.hahsler.net/seriation/reference/Wood.md) :
  Gene Expression Data for Wood Formation in Poplar Trees
- [`Zoo`](https://michael.hahsler.net/seriation/reference/Zoo.md) : Zoo
  Data Set
- [`create_lines_data()`](https://michael.hahsler.net/seriation/reference/create_lines_data.md)
  [`create_ordered_data()`](https://michael.hahsler.net/seriation/reference/create_lines_data.md)
  : Create Simulated Data for Seriation Evaluation
- [`is.robinson()`](https://michael.hahsler.net/seriation/reference/is.robinson.md)
  [`random.robinson()`](https://michael.hahsler.net/seriation/reference/is.robinson.md)
  : Create and Recognize Robinson and Pre-Robinson Matrices

## Helper

- [`LS_swap()`](https://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_insert()`](https://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_reverse()`](https://michael.hahsler.net/seriation/reference/LS.md)
  [`LS_mixed()`](https://michael.hahsler.net/seriation/reference/LS.md)
  : Neighborhood functions for Seriation Method SA
- [`lle()`](https://michael.hahsler.net/seriation/reference/lle.md) :
  Locally Linear Embedding (LLE)
- [`uniscale()`](https://michael.hahsler.net/seriation/reference/uniscale.md)
  [`MDS_stress()`](https://michael.hahsler.net/seriation/reference/uniscale.md)
  [`get_config()`](https://michael.hahsler.net/seriation/reference/uniscale.md)
  [`plot_config()`](https://michael.hahsler.net/seriation/reference/uniscale.md)
  : Fit an Unidimensional Scaling for a Seriation Order
