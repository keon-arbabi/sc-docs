# Basic Workflow

This tutorial walks through a standard single-cell analysis pipeline: feature selection, normalization, dimensionality reduction, clustering, and marker gene identification. It starts with steps from the [Loading and Quality Control](loading_and_qc.md) tutorial using the same Parse Biosciences ~10 million cell PBMC dataset.

## Setup

```python
from single_cell import SingleCell, Timer
import polars as pl
import os, psutil

print(f'{os.cpu_count()} CPUs, '
      f'{psutil.virtual_memory().total / 2**30:.0f} GB RAM')
```
```none
192 CPUs, 755 GB RAM
```

```python
with Timer('Load'):
    sc = SingleCell(
        'Parse_10M_PBMC_cytokines.h5ad', num_threads=-1,
        obs_columns=['sample', 'donor', 'cell_type', 'treatment', 'cytokine'])
```
```none
Load...
Load took 1m 9s
```

```python
with Timer('QC'):
    sc = sc.qc(remove_doublets=True, batch_column='sample',
               subset=False, allow_float=True)
```
```none
QC...
QC took 1s 570ms
```

## Feature selection

{meth}`~single_cell.SingleCell.hvg` selects highly variable genes using Seurat's variance-stabilization approach. It operates on raw counts and must be run before {meth}`~single_cell.SingleCell.normalize`. By default, it selects the top 2,000 genes.

When your data has multiple batches, pass `batch_column` to identify HVGs that are variable across batches rather than just within one:

```python
with Timer('HVG'):
    sc = sc.hvg(batch_column='donor')
```
```none
HVG...
HVG took 2s 280ms
```

This adds a `highly_variable` Boolean column and a `highly_variable_rank` integer column to {attr}`~single_cell.SingleCell.var`. Downstream methods like {meth}`~single_cell.SingleCell.pca` automatically use only the highly variable genes.

## Normalization

{meth}`~single_cell.SingleCell.normalize` log-transforms the counts. The default method, PFlog1pPF ([Booeshaghi et al. 2022](https://biorxiv.org/content/10.1101/2022.05.06.490859v1.full)), applies proportional fitting before and after log1p to correct for composition bias. With `method='logCP10k'`, it matches Seurat's `NormalizeData()`.

```python
with Timer('Normalize'):
    sc = sc.normalize()
```
```none
Normalize...
Normalize took 27s 760ms
```

## PCA

{meth}`~single_cell.SingleCell.pca` computes principal components from the normalized, highly variable genes. The default is 50 PCs, and can be adjusted with `num_PCs`.

```python
with Timer('PCA'):
    sc = sc.pca()
```
```none
PCA...
PCA took 20s 630ms
```

PCs are stored in `obsm['pca']` as a NumPy array of shape `(n_cells, 50)`.

## Nearest neighbors

{meth}`~single_cell.SingleCell.neighbors` computes the 20 nearest neighbors of each cell in PC space, using a fast approximate search based on FAISS's `IndexIVFFlat` strategy. {meth}`~single_cell.SingleCell.shared_neighbors` then builds a shared nearest neighbor (SNN) graph — two cells are connected in proportion to how many neighbors they share. This is the same approach used by Seurat.

```python
with Timer('Neighbors'):
    sc = sc.neighbors().shared_neighbors()
```
```none
Neighbors...
Neighbors took 58s 570ms
```

The nearest-neighbor indices, distances, and SNN graph are stored in `obsm['neighbors']`, `obsm['distances']`, and `obsm['shared_neighbors']`.

:::{note}
If you subset your data after computing neighbors (e.g. via {meth}`~single_cell.SingleCell.filter_obs`), the neighbor graph becomes invalid and must be recomputed. brisc enforces this and will raise an error if you try to use stale neighbors.
:::

## Batch integration

If your dataset contains batch effects, {meth}`~single_cell.SingleCell.harmonize` corrects the PCs using [Harmony](https://github.com/immunogenomics/harmony). Run it after {meth}`~single_cell.SingleCell.pca` and before {meth}`~single_cell.SingleCell.neighbors`:

```python
sc = sc.harmonize(batch_column='donor')
# downstream steps will use obsm['harmony'] instead of obsm['pca']
sc = sc.neighbors(PC_key='harmony').shared_neighbors()
```

brisc's Harmony implementation is parallelized via a nested block strategy and uses a k-means variant called [k-means||](https://arxiv.org/abs/1203.6402) for initialization.

## Clustering

{meth}`~single_cell.SingleCell.cluster` runs Leiden clustering on the SNN graph. The `resolution` parameter controls granularity — higher values produce more clusters. You can pass multiple resolutions to evaluate them in parallel:

```python
with Timer('Cluster'):
    sc = sc.cluster(resolution=[0.25, 0.5, 1.0, 1.5, 2.0])
```
```none
Cluster...
Cluster took 5m 21s
```

Each resolution adds a column to {attr}`~single_cell.SingleCell.obs`: `cluster_0` through `cluster_4` (or a custom name via `cluster_column`). The implementation is based on [GVE-Leiden](https://arxiv.org/abs/2312.13936).

## Embedding

brisc offers three embedding methods for visualization, all of which operate on the PCs and neighbor graph:

- {meth}`~single_cell.SingleCell.pacmap` — [PaCMAP](https://arxiv.org/abs/2012.04456), captures global structure well. Default choice.
- {meth}`~single_cell.SingleCell.localmap` — [LocalMAP](https://arxiv.org/abs/2412.15426), balances local and global structure.
- {meth}`~single_cell.SingleCell.umap` — the standard UMAP algorithm.

```python
with Timer('PaCMAP'):
    sc = sc.pacmap()
```
```none
PaCMAP...
PaCMAP took 1m 8s
```

Embeddings are stored as 2-column NumPy arrays in `obsm` (e.g. `obsm['pacmap']`). Visualize with {meth}`~single_cell.SingleCell.plot_embedding`:

```python
sc.plot_embedding('cell_type', embedding_key='pacmap')
```

## Marker genes

{meth}`~single_cell.SingleCell.find_markers` identifies genes that distinguish each cell type from the rest, using an adaptation of [Fischer and Gillis 2021](https://ncbi.nlm.nih.gov/pmc/articles/PMC8571500). It works on binarized expression (detected vs. not detected), so it gives the same result before or after normalization.

Genes must have a detection rate of at least 25% in the target cell type and at least a 2-fold enrichment relative to every other cell type. Pareto-optimal genes — those not dominated by any other gene on both metrics — are flagged.

```python
with Timer('Markers'):
    markers = sc.find_markers('cell_type')
```
```none
Markers...
Markers took 3s 630ms
```

```python
markers.head()
```
```none
shape: (5, 5)
 cell_type   gene     detection_rate  min_fold_change  pareto
 str         str      f64             f64              bool
 CD4 T       IL7R     0.72            3.1              true
 CD4 T       LTB      0.68            2.8              true
 ...
```

## Pipeline summary

The full pipeline chains together in ~10 minutes on 10 million cells:

```python
sc = (
    SingleCell('data.h5ad', num_threads=-1)
    .qc(remove_doublets=True, batch_column='sample',
        subset=False, allow_float=True)
    .hvg(batch_column='donor')
    .normalize()
    .pca()
    .neighbors()
    .shared_neighbors()
    .cluster(resolution=[0.25, 0.5, 1.0, 1.5, 2.0])
    .pacmap()
)
markers = sc.find_markers('cell_type')
```

| Step | Method | Time | What it does |
|---|---|---|---|
| Feature selection | {meth}`sc.hvg() <single_cell.SingleCell.hvg>` | 2s | Select top 2,000 highly variable genes |
| Normalization | {meth}`sc.normalize() <single_cell.SingleCell.normalize>` | 28s | PFlog1pPF log-normalization |
| PCA | {meth}`sc.pca() <single_cell.SingleCell.pca>` | 21s | 50 principal components |
| Neighbors | {meth}`sc.neighbors() <single_cell.SingleCell.neighbors>` | 59s | 20 nearest neighbors + SNN graph |
| Clustering | {meth}`sc.cluster() <single_cell.SingleCell.cluster>` | 5m 21s | Leiden clustering at multiple resolutions |
| Embedding | {meth}`sc.pacmap() <single_cell.SingleCell.pacmap>` | 1m 8s | 2D PaCMAP embedding |
| Markers | {meth}`sc.find_markers() <single_cell.SingleCell.find_markers>` | 4s | Marker genes per cell type |
