# Integration and Label Transfer

This tutorial covers reference mapping: integrating an annotated reference dataset with a query, then transferring the reference's cell-type labels onto the query. It uses the same Parse Biosciences ~10 million cell PBMC dataset as the other tutorials, treating the PBS-control cells as the annotated reference and the cytokine-treated cells as the query. Because these query cells keep their original labels, we can check the transferred labels against them at the end.

## Loading and quality control

Load the data and run QC as in [Basic Workflow](basic_workflow.md).

```python
from brisc import SingleCell
import polars as pl

sc = SingleCell(
    'Parse_10M_PBMC_cytokines.h5ad',
    obs_columns=['sample', 'donor', 'cell_type', 'treatment', 'cytokine'])\
    .qc(allow_float=True)
```

## Reference and query

{meth}`~brisc.SingleCell.split_by_obs` splits a dataset into one {class}`~brisc.SingleCell` per value of an `obs` column. Reference and query are normally separate datasets; for this self-contained example we split one on `treatment` into the PBS controls (the reference) and cytokine-treated cells (the query).

```python
sc = sc.split_by_obs('treatment')
sc_ref, sc_query = sc['PBS'], sc['cytokine']
print(sc_ref)
print(sc_query)
```
```none
SingleCell dataset in CSR format with 603,928 cells (obs), 40,352 genes (var), and 1,164,409,519 non-zero float32 entries (X)
    obs: _index, sample, donor, cell_type, treatment, cytokine, passed_QC
    var: _index, n_cells
SingleCell dataset in CSR format with 8,839,235 cells (obs), 40,352 genes (var), and 17,313,941,540 non-zero float32 entries (X)
    obs: _index, sample, donor, cell_type, treatment, cytokine, passed_QC
    var: _index, n_cells
```

The reference is far smaller than the query — typical of reference mapping, where a small, carefully annotated dataset labels a much larger one.

## Integration

Integration places the cells of two datasets in one batch-corrected coordinate space, so cells of the same type align regardless of their source dataset. Three steps build that space, each needing both datasets at once: {meth}`~brisc.SingleCell.hvg` picks one shared set of highly variable genes, {meth}`~brisc.SingleCell.pca` builds one shared set of PCs, and {meth}`~brisc.SingleCell.harmonize` removes the batch differences between them with [Harmony](https://github.com/immunogenomics/harmony), storing the result in `obsm['harmony']`. So each takes the other dataset as an argument and returns both.

Normalization is the exception — it treats each cell independently, so the datasets are normalized separately.

```python
sc_ref, sc_query = sc_ref.hvg(sc_query, batch_column='donor')
sc_ref = sc_ref.normalize()
sc_query = sc_query.normalize()
sc_ref, sc_query = sc_ref.pca(sc_query)
sc_ref, sc_query = sc_ref.harmonize(sc_query)
```

```none
40,352 genes are present in every dataset.
Initialization is complete: objective = 3932457.75
Completed 1 of 10 iterations: objective = 2163470.00 (k-means error = 2726398.25, entropy term = -1872335.38, diversity penalty = 1309407.12)
Completed 2 of 10 iterations: objective = 2146725.00 (k-means error = 2740209.50, entropy term = -1902801.25, diversity penalty = 1309316.88)
Reached convergence after 2 iterations
```

## Label transfer

{meth}`~brisc.SingleCell.label_transfer_from` transfers cell-type labels from the reference to the query. For each query cell, it finds the `num_neighbors` (default 20) nearest reference cells in the shared Harmony embedding and assigns the most common reference label; the fraction of those neighbors that agree becomes a confidence score.

```python
sc_query = sc_query.label_transfer_from(
    sc_ref, 'cell_type', cell_type_column='cell_type_transferred')
```

This adds `cell_type_transferred` and `cell_type_transferred_confidence` to {attr}`~brisc.SingleCell.obs`:

```python
print(sc_query.obs.select('cell_type', 'cell_type_transferred',
                          'cell_type_transferred_confidence').head(10))
```
```none
shape: (10, 3)
┌────────────┬───────────────────────┬──────────────────────────────────┐
│ cell_type  ┆ cell_type_transferred ┆ cell_type_transferred_confidence │
│ ---        ┆ ---                   ┆ ---                              │
│ enum       ┆ enum                  ┆ f32                              │
╞════════════╪═══════════════════════╪══════════════════════════════════╡
│ CD8 Naive  ┆ CD8 Naive             ┆ 1.0                              │
│ B Naive    ┆ B Naive               ┆ 0.85                             │
│ CD14 Mono  ┆ CD14 Mono             ┆ 1.0                              │
│ CD14 Mono  ┆ CD14 Mono             ┆ 0.9                              │
│ CD4 Naive  ┆ CD4 Naive             ┆ 0.95                             │
│ CD8 Naive  ┆ CD4 Memory            ┆ 0.45                             │
│ NK         ┆ NK                    ┆ 1.0                              │
│ CD4 Memory ┆ CD4 Memory            ┆ 0.95                             │
│ NK         ┆ NK                    ┆ 1.0                              │
│ cDC        ┆ cDC                   ┆ 1.0                              │
└────────────┴───────────────────────┴──────────────────────────────────┘
```

Most calls are confident; the low-confidence row (0.45) is a CD8 Naive cell labeled CD4 Memory, and the confidence score is what flags these uncertain transfers so you can filter on them.

:::{dropdown} Runner-up labels
Pass `next_best=True` to also record each cell's runner-up label and its confidence, in `next_best_cell_type_transferred` and `next_best_cell_type_transferred_confidence`. This helps when a cell sits between two similar types:

```python
sc_query = sc_query.label_transfer_from(
    sc_ref, 'cell_type', cell_type_column='cell_type_transferred', next_best=True)
print(sc_query.obs.select(
    'cell_type', 'cell_type_transferred', 'cell_type_transferred_confidence',
    'next_best_cell_type_transferred',
    'next_best_cell_type_transferred_confidence').head(10))
```
```none
shape: (10, 5)
┌────────────┬───────────────────────┬───────────────────────────┬───────────────────────────┬───────────────────────────┐
│ cell_type  ┆ cell_type_transferred ┆ cell_type_transferred_con ┆ next_best_cell_type_trans ┆ next_best_cell_type_trans │
│ ---        ┆ ---                   ┆ fiden…                    ┆ ferre…                    ┆ ferre…                    │
│ enum       ┆ enum                  ┆ ---                       ┆ ---                       ┆ ---                       │
│            ┆                       ┆ f32                       ┆ enum                      ┆ f32                       │
╞════════════╪═══════════════════════╪═══════════════════════════╪═══════════════════════════╪═══════════════════════════╡
│ CD8 Naive  ┆ CD8 Naive             ┆ 1.0                       ┆ CD4 Memory                ┆ 0.0                       │
│ B Naive    ┆ B Naive               ┆ 0.85                      ┆ B Intermediate/Memory     ┆ 0.15                      │
│ CD14 Mono  ┆ CD14 Mono             ┆ 1.0                       ┆ CD4 Memory                ┆ 0.0                       │
│ CD14 Mono  ┆ CD14 Mono             ┆ 0.9                       ┆ CD16 Mono                 ┆ 0.05                      │
│ CD4 Naive  ┆ CD4 Naive             ┆ 0.95                      ┆ CD4 Memory                ┆ 0.05                      │
│ CD8 Naive  ┆ CD4 Memory            ┆ 0.45                      ┆ CD8 Naive                 ┆ 0.45                      │
│ NK         ┆ NK                    ┆ 1.0                       ┆ CD4 Memory                ┆ 0.0                       │
│ CD4 Memory ┆ CD4 Memory            ┆ 0.95                      ┆ MAIT                      ┆ 0.05                      │
│ NK         ┆ NK                    ┆ 1.0                       ┆ CD4 Memory                ┆ 0.0                       │
│ cDC        ┆ cDC                   ┆ 1.0                       ┆ CD4 Memory                ┆ 0.0                       │
└────────────┴───────────────────────┴───────────────────────────┴───────────────────────────┴───────────────────────────┘
```

Here the 0.45 cell splits evenly between CD4 Memory and its runner-up, CD8 Naive — which is its true label.
:::

## Validation

Because the query carries ground-truth labels, we can measure how well the transfer recovered them.

```python
match = (pl.col('cell_type').cast(pl.String) ==
         pl.col('cell_type_transferred').cast(pl.String))
print(f'overall accuracy: {sc_query.obs.select(match.mean()).item():.1%}')
print(sc_query.obs.group_by('cell_type').agg(
    n_cells=pl.len(),
    accuracy=match.mean(),
    mean_confidence=pl.col('cell_type_transferred_confidence').mean())
    .sort('accuracy', descending=True))
```
```none
overall accuracy: 88.1%
shape: (18, 4)
┌───────────────────────┬─────────┬──────────┬─────────────────┐
│ cell_type             ┆ n_cells ┆ accuracy ┆ mean_confidence │
│ ---                   ┆ ---     ┆ ---      ┆ ---             │
│ enum                  ┆ u32     ┆ f64      ┆ f32             │
╞═══════════════════════╪═════════╪══════════╪═════════════════╡
│ CD14 Mono             ┆ 1443470 ┆ 0.993401 ┆ 0.988664        │
│ Plasmablast           ┆ 3908    ┆ 0.973132 ┆ 0.969703        │
│ cDC                   ┆ 102056  ┆ 0.970987 ┆ 0.974583        │
│ NK                    ┆ 475494  ┆ 0.970096 ┆ 0.956886        │
│ pDC                   ┆ 17561   ┆ 0.966346 ┆ 0.97694         │
│ B Naive               ┆ 547712  ┆ 0.965303 ┆ 0.959042        │
│ HSPC                  ┆ 15324   ┆ 0.938593 ┆ 0.967753        │
│ CD16 Mono             ┆ 212214  ┆ 0.935381 ┆ 0.945373        │
│ B Intermediate/Memory ┆ 279140  ┆ 0.907025 ┆ 0.920662        │
│ CD4 Naive             ┆ 1603760 ┆ 0.889186 ┆ 0.849431        │
│ ILC                   ┆ 7780    ┆ 0.879692 ┆ 0.942018        │
│ CD8 Naive             ┆ 573696  ┆ 0.879201 ┆ 0.865165        │
│ NK CD56bright         ┆ 116694  ┆ 0.856462 ┆ 0.918812        │
│ CD8 Memory            ┆ 684115  ┆ 0.846543 ┆ 0.821517        │
│ CD4 Memory            ┆ 2164338 ┆ 0.806086 ┆ 0.818958        │
│ MAIT                  ┆ 289824  ┆ 0.792943 ┆ 0.840404        │
│ Treg                  ┆ 156180  ┆ 0.791087 ┆ 0.864583        │
│ NKT                   ┆ 145969  ┆ 0.433667 ┆ 0.727498        │
└───────────────────────┴─────────┴──────────┴─────────────────┘
```

Common, distinct types transfer almost perfectly (CD14 Mono 99%, NK 97%, B Naive 97%), while rare or closely related types are harder — NKT (43%) is mostly absorbed into the neighboring NK and CD8 populations, and its low mean confidence (0.73) reflects that.

## Pipeline summary

The full reference-mapping pipeline:

```python
sc = SingleCell('data.h5ad').qc(allow_float=True)
sc = sc.split_by_obs('treatment')
sc_ref, sc_query = sc['PBS'], sc['cytokine']
sc_ref, sc_query = sc_ref.hvg(sc_query, batch_column='donor')
sc_ref = sc_ref.normalize()
sc_query = sc_query.normalize()
sc_ref, sc_query = sc_ref.pca(sc_query)
sc_ref, sc_query = sc_ref.harmonize(sc_query)
sc_query = sc_query.label_transfer_from(
    sc_ref, 'cell_type', cell_type_column='cell_type_transferred')
```

| Step | Method | What it does |
|---|---|---|
| Split | {meth}`sc.split_by_obs('treatment') <brisc.SingleCell.split_by_obs>` | Separate the reference (PBS) and query (cytokine) |
| Feature selection | {meth}`sc_ref.hvg(sc_query) <brisc.SingleCell.hvg>` | Highly variable genes across both datasets |
| Normalization | {meth}`sc_ref.normalize() <brisc.SingleCell.normalize>` | log1pPF log-normalization, per dataset |
| PCA | {meth}`sc_ref.pca(sc_query) <brisc.SingleCell.pca>` | A shared principal-component space |
| Integration | {meth}`sc_ref.harmonize(sc_query) <brisc.SingleCell.harmonize>` | Harmony-corrected embedding in `obsm['harmony']` |
| Label transfer | {meth}`sc_query.label_transfer_from(sc_ref) <brisc.SingleCell.label_transfer_from>` | Transfer labels via nearest neighbors in Harmony space |
