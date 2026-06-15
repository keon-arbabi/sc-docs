# Essential Commands

A cheatsheet of the commands for inspecting, subsetting, and reshaping a dataset.

## Inspecting

```python
# summary: dimensions, X dtype, and all obs/var/obsm/… keys
print(sc)

# shape (n_cells, n_genes) / cell count
sc.shape
len(sc)

# column names / columns with their dtypes
sc.obs.columns
sc.obs.schema

# first cell / first gene, one column per line
sc.peek_obs()
sc.peek_var()
```

Inspect or partially load a file without reading the matrix:

```python
# a file's dimensions and obs/var columns
SingleCell.ls('data.h5ad')

# read individual slots without loading X
obs = SingleCell.read_obs('data.h5ad')
var = SingleCell.read_var('data.h5ad')
obsm = SingleCell.read_obsm('data.h5ad')
varm = SingleCell.read_varm('data.h5ad')
uns = SingleCell.read_uns('data.h5ad')
```

## Accessing the pieces

```python
# cell / gene tables (polars DataFrames)
sc.obs
sc.var

# one column (a polars Series)
sc.obs['cell_type']

# cell barcodes / gene names (first columns of obs / var)
sc.obs_names
sc.var_names

# the count matrix (sparse, cells × genes, CSR)
sc.X

# one gene across cells / one cell across genes (dense 1D array)
sc.gene('MALAT1')
sc.cell(sc.obs_names[0])

# pipeline outputs (after PCA / clustering — see Basic Workflow)
sc.obsm['pca']                # PCA embedding (NumPy array)
sc.obsp['shared_neighbors']   # cell–cell graph (sparse)
sc.uns['QCed']                # pipeline flags
```

## Subsetting cells and genes

By metadata — `X` stays in sync:

```python
# one condition / the equality shorthand
sc = sc.filter_obs(pl.col('cell_type') == 'CD14 Mono')
sc = sc.filter_obs(cell_type='CD14 Mono')

# membership in a set
sc = sc.filter_obs(pl.col('cytokine').is_in(['IFN-gamma', 'PBS']))

# several predicates (combined with AND)
sc = sc.filter_obs(pl.col('treatment') == 'PBS', pl.col('cell_type') == 'CD14 Mono')

# genes instead of cells
sc = sc.filter_var(pl.col('n_cells') > 100)
```

By position or name:

```python
# by position / by barcode (all genes)
sc[0]
sc[sc.obs_names[0]]

# one gene / several genes (all cells)
sc[:, 'MALAT1']
sc[:, ['CD3D', 'CD8A', 'MS4A1']]

# first 1000 cells / a Boolean mask
sc[:1000]
sc[sc.obs['cell_type'] == 'CD14 Mono']
```

Split or sample (`split_by_var` / `subsample_var` do the same for genes):

```python
# split into a dict of datasets, keyed by column value
parts = sc.split_by_obs('treatment')

# 10k random cells / 10% of each cell type
sc = sc.subsample_obs(n=10_000)
sc = sc.subsample_obs(fraction=0.1, by_column='cell_type')
```

## Editing metadata

```python
# add or overwrite a column from an expression
sc = sc.with_columns_obs(
    stimulated=pl.col('treatment') == 'cytokine')

# add a multi-branch column
sc = sc.with_columns_obs(
    lineage=pl.when(pl.col('cell_type').cast(pl.String).str.contains('Mono|DC'))
              .then(pl.lit('myeloid'))
              .otherwise(pl.lit('lymphoid')))

# rename / drop / keep / retype columns
sc = sc.rename_obs({'treatment': 'group'})
sc = sc.drop_obs('cytokine')
sc = sc.select_obs('cell_type', 'sample')   # keeps obs_names automatically
sc = sc.cast_obs({'donor': pl.String})
```

Left-join an external table onto `obs` (per-donor or per-sample annotations). The join keys must share a dtype:

```python
donor_meta = pl.DataFrame({'donor': ['Donor1', 'Donor2'], 'age': [29, 34]})

# sc.obs['donor'] is an Enum, so match its dtype before joining
donor_meta = donor_meta.cast({'donor': sc.obs.schema['donor']})
sc = sc.join_obs(donor_meta, on='donor')

# the gene-side twin
sc = sc.join_var(gene_meta, on='gene_id')
```

Set or de-duplicate the names:

```python
# set a column as the cell names (String, Enum, Categorical, or integer)
sc = sc.set_obs_names('sample')

# de-duplicate repeated cell / gene names (appends -1, -2, …)
sc = sc.make_obs_names_unique()
sc = sc.make_var_names_unique()
```

## The count matrix

```python
# the count matrix (sparse, CSR — fast for per-cell access)
sc.X

# storage format: 'csr' or 'csc'
sc.X.format

# convert layout: CSC for per-gene access, CSR for cell-wise work
sc = sc.tocsc()
sc = sc.tocsr()

# change the matrix dtype / drop X to save memory
sc = sc.cast_X('float32')
sc = sc.drop_X()

# set the default thread count (-1 = all cores)
sc = sc.set_num_threads(8)
```

## Summarizing

```python
# cells per type / per sample
sc.obs['cell_type'].value_counts(sort=True)
sc.obs.group_by('sample').agg(n=pl.len())

# several stats per cell type at once
sc.obs.group_by('cell_type').agg(
    n=pl.len(),
    n_donors=pl.col('donor').n_unique())

# columns constant within each sample
sc.get_sample_covariates(ID_column='sample')

# add num_counts / num_genes / mito_fraction, then summarize
sc = sc.qc_metrics(allow_float=True)
sc.obs.group_by('cell_type').agg(pl.col('num_counts').median())

# cross-tabulate two categoricals as a heatmap
sc.plot_heatmap('cell_type', 'donor', 'composition.png')
```

## Combining and saving

```python
from brisc import concat_obs, concat_var

# stack cells (same genes), optionally labeling each cell's source
combined = concat_obs([sc_a, sc_b])
combined = concat_obs([sc_a, sc_b], dataset_column='batch')

# intersect genes and columns first, for mismatched datasets
combined = concat_obs([sc_a, sc_b], flexible=True)

# stack genes (same cells, e.g. RNA + protein)
combined = concat_var([sc_rna, sc_adt])

# shallow copy (default) — shares X and the obsm/obsp/uns arrays; fine unless you mutate them in place
sc_copy = sc.copy()

# deep copy — also duplicates X and those arrays, for a fully independent dataset
sc_deep = sc.copy(deep=True)

# save to any format (see Interoperability)
sc.save('output.h5ad')
```

## Pseudobulk

{meth}`~brisc.SingleCell.pseudobulk` returns a {class}`~brisc.Pseudobulk` — three slots (`X`, `obs`, `var`), each a dict keyed by cell type — with the same `*_obs` / `*_var` methods:

```python
# index by cell type (a single-cell-type Pseudobulk)
pb['CD14 Mono']

# the slot dicts, each keyed by cell type
pb.X
pb.obs
pb.var

# one cell type's matrix / sample table / gene table
pb.X['CD14 Mono']
pb.obs['CD14 Mono']
pb.var['CD14 Mono']

# sample IDs / gene names / shapes, per cell type
pb.obs_names
pb.var_names
pb.shape

# the cell types / their (X, obs, var) triples
pb.keys()
pb.values()

# iterate over the cell types
for cell_type in pb:
    ...

# iterate (cell_type, (X, obs, var)) pairs
for cell_type, (X, obs, var) in pb.items():
    ...

# iterate one slot at a time, across cell types
for obs in pb.iter_obs():
    ...
for var in pb.iter_var():
    ...

# any *_obs method, applied within every cell type
pb = pb.filter_obs(pl.col('treatment') == 'PBS')

# operate on the cell types themselves
pb = pb.select_cell_types(['CD14 Mono', 'NK'])
pb = pb.rename_cell_types({'CD14 Mono': 'Mono'})

# one wide table: row per (sample, cell type), gene columns
pb.to_df()

# normalize the expression, then regress out covariates
pb = pb.library_size()
pb = pb.log_CPM()
pb = pb.regress_out('~ donor + log2(num_cells) + log2(library_size)')
```

This `log_CPM` / `regress_out` path is separate from {meth}`~brisc.Pseudobulk.DE`, which needs raw counts and so can't run after `log_CPM`.

## DE results

{meth}`~brisc.Pseudobulk.DE` returns a {class}`~brisc.DE` whose `table` is a polars DataFrame you filter and sort normally:

```python
# one row per gene per cell type
de.table

# significant genes (FDR < 0.05) / hit counts per cell type
de.get_hits(significance_column='FDR', threshold=0.05)
de.get_num_hits()

# filter and sort like any polars DataFrame
de.table.filter(pl.col('cell_type') == 'CD14 Mono').sort('FDR').head(20)

# save, then reload
de.save('de_results')
de = DE('de_results')
```

## Custom operations

For anything not wrapped by a brisc method, reach for `pipe`:

```python
# any polars op on obs (or var), keeping X in sync
sc = sc.pipe_obs(lambda df: df.sort('cell_type'))

# transform X / the whole dataset
sc = sc.pipe_X(lambda X: X.sqrt())
sc = sc.pipe(my_function)
```

A {class}`~brisc.Pseudobulk` has the same hooks; `map_obs` / `map_X` apply the function within each cell type, while `pipe` takes the whole dataset:

```python
pb = pb.map_obs(lambda df: df.sort('sample'))
pb = pb.map_X(lambda X: np.log1p(X))
pb = pb.pipe(my_function)
```
