# Data Manipulation

This tutorial is a cheatsheet for inspecting, subsetting, and reshaping a dataset — the polars equivalents of the metadata operations familiar from scanpy and Seurat. brisc keeps cell and gene metadata in [polars](https://pola.rs) DataFrames rather than pandas, so the tables below map the concepts across. It uses the same Parse Biosciences ~10 million cell PBMC dataset as the other tutorials; snippets use the loaded dataset `sc` (see [Basic Workflow](basic_workflow.md)).

## The SingleCell object

A {class}`~brisc.SingleCell` object bundles a count matrix with its cell and gene metadata:

- {attr}`~brisc.SingleCell.X` — the counts, a SciPy sparse array of shape `(cells, genes)`, in CSR format.
- {attr}`~brisc.SingleCell.obs` — one metadata row per **cell**, a polars DataFrame.
- {attr}`~brisc.SingleCell.var` — one metadata row per **gene**, a polars DataFrame.
- {attr}`~brisc.SingleCell.obsm` / {attr}`~brisc.SingleCell.obsp` — per-cell matrices (embeddings like `obsm['pca']`) and cell–cell graphs (`obsp['shared_neighbors']`).
- {attr}`~brisc.SingleCell.uns` — unstructured odds and ends (flags, parameters).

Two important things to remember:

- **No row index.** polars DataFrames have no index. Cells and genes are identified by the **first column** of `obs` and `var`, surfaced as {attr}`~brisc.SingleCell.obs_names` and {attr}`~brisc.SingleCell.var_names`.
- **Immutable.** Every method returns a *new* `SingleCell`; nothing changes in place. Assign the result to keep it — `sc = sc.filter_obs(...)`.

If you are coming from AnnData or Seurat:

| Concept | AnnData | Seurat | brisc |
|---|---|---|---|
| Counts | `adata.X` | `obj[["RNA"]]$counts` | `sc.X` (CSR, cells × genes) |
| Cell metadata | `adata.obs` (pandas) | `obj@meta.data` | `sc.obs` (polars) |
| Gene metadata | `adata.var` (pandas) | feature metadata | `sc.var` (polars) |
| Cell names | `adata.obs_names` | `colnames(obj)` | `sc.obs_names` (first `obs` column) |
| Gene names | `adata.var_names` | `rownames(obj)` | `sc.var_names` (first `var` column) |
| Embeddings | `adata.obsm['X_pca']` | `obj@reductions` | `sc.obsm['pca']` |
| Cell–cell graphs | `adata.obsp` | `obj@graphs` | `sc.obsp` |
| Unstructured | `adata.uns` | `obj@misc` | `sc.uns` |

## Working with obs and var

`sc.obs` and `sc.var` are ordinary polars DataFrames. Read them directly for anything that does not change the dataset:

```python
# the whole cell table
sc.obs
# one column → a polars Series
sc.obs['cell_type']
# count occurrences of each value, most frequent first
sc.obs['cell_type'].value_counts(sort=True)
# number of distinct donors
sc.obs['donor'].n_unique()
```

To change the dataset and keep `X` aligned, use the matching `*_obs` / `*_var` methods. Each wraps the polars verb of the same name — `filter`, `select`, `with_columns`, `rename`, `drop`, `join` — applied to the table and the count matrix together, so {meth}`sc.filter_obs(...) <brisc.SingleCell.filter_obs>` is `sc.obs.filter(...)` with `X` subset to match. Use `sc.obs` directly when you want a table back, a `*_obs` method when you want a dataset.

Coming from pandas, the everyday translations are:

| Task | pandas / AnnData | brisc / polars |
|---|---|---|
| A column | `adata.obs['cell_type']` | `sc.obs['cell_type']` |
| Count categories | `adata.obs['cell_type'].value_counts()` | `sc.obs['cell_type'].value_counts()` |
| Filter cells | `adata[adata.obs['cell_type'] == 'CD14 Mono']` | `sc.filter_obs(pl.col('cell_type') == 'CD14 Mono')` |
| Add a column | `adata.obs['stim'] = adata.obs['treatment'] == 'cytokine'` | `sc = sc.with_columns_obs(stim=pl.col('treatment') == 'cytokine')` |
| Rename | `adata.obs.rename(columns={'treatment': 'group'})` | `sc = sc.rename_obs({'treatment': 'group'})` |
| Group + mean | `adata.obs.groupby('sample')['num_counts'].mean()` | `sc.obs.group_by('sample').agg(pl.col('num_counts').mean())` |
| Genes by name | `adata[:, ['CD3D', 'CD8A']]` | `sc[:, ['CD3D', 'CD8A']]` |

Three things to unlearn from pandas:

- **Expressions, not masks.** Conditions are built with `pl.col` — `pl.col('num_counts') > 1000` — and combined with `&`, `|`, `~`, each side parenthesized, not `and`/`or`.
- **Assign the result.** polars and brisc are both immutable; `sc.with_columns_obs(...)` is a no-op unless you write `sc = sc.with_columns_obs(...)`.
- **No `.loc`/`.iloc`.** There is no index to align on — pick columns by name, and rows by filtering or positional slicing (`sc.obs[5:10]`).

## Inspecting

```python
# dimensions, X dtype, and every obs/var/obsm/... key
print(sc)
# (n_cells, n_genes)
sc.shape
# n_cells
len(sc)
# cell-metadata column names
sc.obs.columns
# columns with their polars dtypes
sc.obs.schema
# first cell, one column per line (no truncation)
sc.peek_obs()
# first gene
sc.peek_var()
```

Inspect a file's contents — dimensions and available `obs`/`var` columns — without reading the matrix:

```python
SingleCell.ls('Parse_10M_PBMC_cytokines.h5ad')
```

Or load only the metadata (no `X`):

```python
obs = SingleCell.read_obs('Parse_10M_PBMC_cytokines.h5ad')
var = SingleCell.read_var('Parse_10M_PBMC_cytokines.h5ad')
```

## Accessing the pieces

```python
# cell and gene tables (polars DataFrames)
sc.obs
sc.var
# one column → a polars Series
sc.obs['cell_type']
# cell barcodes and gene names (first columns of obs / var)
sc.obs_names
sc.var_names
# the counts: a sparse array, cells × genes (CSR)
sc.X
# one gene across all cells → dense 1D array
sc.gene('MALAT1')
# one cell across all genes → dense 1D array
sc.cell(sc.obs_names[0])
```

After you run PCA, clustering, and embeddings (see [Basic Workflow](basic_workflow.md)), their outputs live in `obsm`, `obsp`, and `uns`:

```python
# available cell embeddings
sc.obsm.keys()
# PCA coordinates (NumPy array)
sc.obsm['pca']
# the shared-nearest-neighbor graph (sparse)
sc.obsp['shared_neighbors']
# flags set by the pipeline
sc.uns['QCed']
```

## Subsetting cells and genes

By metadata — `X` stays in sync:

```python
# one condition
sc.filter_obs(pl.col('cell_type') == 'CD14 Mono')
# shorthand for equality
sc.filter_obs(cell_type='CD14 Mono')
# membership in a set
sc.filter_obs(pl.col('cytokine').is_in(['IFN-gamma', 'PBS']))
# several predicates combine with AND
sc.filter_obs(pl.col('treatment') == 'PBS', pl.col('cell_type') == 'CD14 Mono')
# filter genes instead
sc.filter_var(pl.col('n_cells') > 100)
```

By position or name:

```python
# first cell (all genes)
sc[0]
# a cell by barcode
sc[sc.obs_names[0]]
# one gene (all cells)
sc[:, 'MALAT1']
# several genes
sc[:, ['CD3D', 'CD8A', 'MS4A1']]
# first 1000 cells
sc[:1000]
# a Boolean mask also works
sc[sc.obs['cell_type'] == 'CD14 Mono']
```

Split into a dictionary of datasets, or take a random sample:

```python
# split into {'PBS': SingleCell, 'cytokine': SingleCell}
by_treatment = sc.split_by_obs('treatment')
# 10k random cells
sc.subsample_obs(n=10_000)
# 10% of each cell type
sc.subsample_obs(fraction=0.1, by_column='cell_type')
```

## Editing metadata

Add or overwrite columns with `pl` expressions:

```python
# a Boolean column from a comparison
sc = sc.with_columns_obs(stimulated=pl.col('treatment') == 'cytokine')
# a categorical column from a multi-branch condition
sc = sc.with_columns_obs(
    lineage=pl.when(pl.col('cell_type').cast(pl.String).str.contains('Mono|DC'))
              .then(pl.lit('myeloid'))
              .otherwise(pl.lit('lymphoid')))
```

Rename, drop, keep, or retype columns:

```python
# rename a column
sc = sc.rename_obs({'treatment': 'group'})
# drop a column
sc = sc.drop_obs('cytokine')
# keep only these columns (obs_names is kept automatically)
sc = sc.select_obs('cell_type', 'sample')
# change a column's dtype
sc = sc.cast_obs({'donor': pl.String})
```

Left-join an external table onto `obs` — per-donor or per-sample annotations — with the same logic as polars' `DataFrame.join`. The join keys must share a dtype, so cast the external column to match the one in `obs`:

```python
donor_meta = pl.DataFrame({'donor': ['Donor1', 'Donor2'], 'age': [29, 34]})
# obs['donor'] is an Enum, so match its dtype before joining
donor_meta = donor_meta.cast({'donor': sc.obs.schema['donor']})
sc = sc.join_obs(donor_meta, on='donor')
```

{meth}`~brisc.SingleCell.join_var` does the same on the gene side, e.g. attaching biotype or pathway annotations.

Cell and gene names are just the first columns of `obs` and `var`. Set a different column as the names with {meth}`~brisc.SingleCell.set_obs_names` (and {meth}`~brisc.SingleCell.set_var_names` for genes), or make repeated names unique — appending `-1` to the second occurrence, `-2` to the third, and so on — with {meth}`~brisc.SingleCell.make_obs_names_unique` / {meth}`~brisc.SingleCell.make_var_names_unique`:

```python
# promote a column to be the cell names (String, Enum, Categorical, or integer)
sc = sc.set_obs_names('sample')
# disambiguate duplicates: 'CD8', 'CD8' → 'CD8', 'CD8-1'
sc = sc.make_var_names_unique()
```

## Summarizing and grouping

Summaries read straight off `sc.obs` (no dataset returned):

```python
# cells per type
sc.obs['cell_type'].value_counts(sort=True)
# cells per sample
sc.obs.group_by('sample').agg(n=pl.len())
# several stats per cell type at once
sc.obs.group_by('cell_type').agg(
    n=pl.len(),
    n_donors=pl.col('donor').n_unique())
# columns that are constant within each sample
sc.get_sample_covariates(ID_column='sample')
```

{meth}`~brisc.SingleCell.qc_metrics` adds `num_counts`, `num_genes`, and `mito_fraction` to `obs`, which you can then summarize like any column:

```python
sc = sc.qc_metrics(allow_float=True)
sc.obs.group_by('cell_type').agg(pl.col('num_counts').median())
```

Cross-tabulate two categorical columns as a heatmap:

```python
sc.plot_heatmap('cell_type', 'donor', 'composition.png')
```

## The count matrix

```python
# the sparse counts — CSR (the default) is fast for per-cell access
sc.X
# the storage layout: 'csr' or 'csc'
sc.X.format
# convert to CSC for fast per-gene access...
sc = sc.tocsc()
# ...and back to CSR, the better layout for cell-wise work like pseudobulking
sc = sc.tocsr()
# change the matrix dtype
sc = sc.cast_X('float32')
# one gene's counts across cells → dense vector
sc.gene('CD3D')
```

Per-cell totals come from {meth}`~brisc.SingleCell.qc_metrics` (`num_counts`), or straight from `X`:

```python
import numpy as np
# total counts per cell
np.asarray(sc.X.sum(axis=1)).ravel()
```

Normalizing and log-transforming counts is part of the analysis workflow — see {meth}`~brisc.SingleCell.normalize` in [Basic Workflow](basic_workflow.md).

Most operations run multithreaded. Change a dataset's default thread count with {meth}`~brisc.SingleCell.set_num_threads` (`-1` uses all cores), or assign the `num_threads` property in place:

```python
# 8 threads for this dataset's subsequent operations
sc = sc.set_num_threads(8)
```

## Combining datasets

Stack datasets cell-wise or gene-wise:

```python
from brisc import concat_obs, concat_var

# stack cells (same genes)
concat_obs([sc_a, sc_b])
# ...and label each cell's source dataset
concat_obs([sc_a, sc_b], dataset_column='batch')
# intersect genes and columns first, for mismatched datasets
concat_obs([sc_a, sc_b], flexible=True)
# stack genes (same cells, e.g. RNA + protein)
concat_var([sc_rna, sc_adt])
```

{meth}`~brisc.SingleCell.copy` duplicates a dataset. `obs` and `var` are always shared (polars DataFrames are immutable, so this is safe); by default `X` and the NumPy arrays in `obsm`/`obsp`/`uns` are shared too, so in-place edits to them would affect both copies. Pass `deep=True` to copy those arrays as well:

```python
# shallow — shares X and the underlying arrays
sc2 = sc.copy()
# fully independent copy
sc2 = sc.copy(deep=True)
```

## Pseudobulk objects

{meth}`~brisc.SingleCell.pseudobulk` returns a {class}`~brisc.Pseudobulk` — one matrix plus `obs` and `var` per cell type. It behaves like a dictionary keyed by cell type and mirrors the same `*_obs` manipulation methods. Given a pseudobulk (see [Differential Expression](differential_expression.md)):

```python
# index by cell type → a single-cell-type Pseudobulk
pb['CD14 Mono']
# that cell type's sample table (polars DataFrame)
pb.obs['CD14 Mono']
# the cell types
pb.keys()
# iterate over cell types
for cell_type, (X, obs, var) in pb.items():
    ...
# a filter applied within every cell type
pb.filter_obs(pl.col('treatment') == 'PBS')
# one wide table: row per (sample, cell type), gene columns
pb.to_df()
```

## DE results

{meth}`~brisc.Pseudobulk.DE` returns a {class}`~brisc.DE` whose {attr}`~brisc.DE.table` is a polars DataFrame — filter and sort it like any other:

```python
# one row per gene per cell type
de.table
# significant genes (FDR < 0.05)
de.get_hits(significance_column='FDR', threshold=0.05)
# number of hits per cell type
de.get_num_hits()
# filter and sort the table like any polars DataFrame
de.table.filter(pl.col('cell_type') == 'CD14 Mono').sort('FDR').head(20)
```

## Escape hatches

For anything not wrapped by a brisc method, use `pipe`:

```python
# run any polars op on obs, keeping X in sync
sc.pipe_obs(lambda df: df.sort('cell_type'))
# run any function on the whole dataset
sc.pipe(my_function)
```

## Summary

| Operation | Method |
|---|---|
| Inspect a file without loading | {meth}`SingleCell.ls(file) <brisc.SingleCell.ls>` |
| Load only metadata | {meth}`SingleCell.read_obs(file) <brisc.SingleCell.read_obs>` |
| Dimensions | {attr}`sc.shape <brisc.SingleCell.shape>`, `len(sc)` |
| A metadata column | `sc.obs['cell_type']` |
| Count categories | `sc.obs['cell_type'].value_counts()` |
| Filter cells | {meth}`sc.filter_obs() <brisc.SingleCell.filter_obs>` |
| Filter genes | {meth}`sc.filter_var() <brisc.SingleCell.filter_var>` |
| Subset by name or position | `sc[cells, genes]` |
| Split into groups | {meth}`sc.split_by_obs() <brisc.SingleCell.split_by_obs>` |
| Subsample | {meth}`sc.subsample_obs() <brisc.SingleCell.subsample_obs>` |
| Add or edit a column | {meth}`sc.with_columns_obs() <brisc.SingleCell.with_columns_obs>` |
| Rename / drop / keep columns | {meth}`sc.rename_obs() <brisc.SingleCell.rename_obs>`, {meth}`drop_obs() <brisc.SingleCell.drop_obs>`, {meth}`select_obs() <brisc.SingleCell.select_obs>` |
| Merge external metadata | {meth}`sc.join_obs() <brisc.SingleCell.join_obs>` |
| One gene's counts | {meth}`sc.gene() <brisc.SingleCell.gene>` |
| Convert matrix format | {meth}`sc.tocsc() <brisc.SingleCell.tocsc>` |
| Group and summarize | `sc.obs.group_by(...).agg(...)` |
| Combine datasets | {func}`concat_obs() <brisc.concat_obs>` |
| Aggregate to pseudobulk | {meth}`sc.pseudobulk() <brisc.SingleCell.pseudobulk>` |
| Arbitrary operation | {meth}`sc.pipe_obs() <brisc.SingleCell.pipe_obs>` |
