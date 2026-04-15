# Loading and Quality Control

## Dataset

This tutorial uses a ~10 million cell PBMC cytokine stimulation dataset from [Parse Biosciences](https://www.parsebiosciences.com/datasets). Cryopreserved PBMCs from twelve healthy donors were seeded at 1 million cells per well in 96-well plates -- one plate per donor -- and treated with 90 different cytokines or PBS control for 24 hours, yielding 1,092 experimental conditions. Fixed samples were processed in the Parse GigaLab and sequenced on Ultima Genomics at ~31,000 mean reads per cell. The dataset is already processed, but we will walk through quality-control filtering and doublet detection nonetheless.

## Setup

```python
from single_cell import SingleCell, Timer
from subprocess import run
import os, psutil

print(f'{os.cpu_count()} CPUs, '
      f'{psutil.virtual_memory().total / 2**30:.0f} GB RAM')
```
```none
192 CPUs, 755 GB RAM
```

```python
run('wget -nc https://parse-wget.s3.us-west-2.amazonaws.com/10m/'
    'Parse_10M_PBMC_cytokines.h5ad',
    shell=True)
```

## Inspecting the file

Before loading, you can inspect the `.h5ad` file's structure with {meth}`~single_cell.SingleCell.ls` to see what's inside. Useful if you want to load only selected metadata columns to reduce loading time and memory usage.

```python
SingleCell.ls('Parse_10M_PBMC_cytokines.h5ad')
```
```none
X: 9,697,974 × 40,352 sparse array with 18,830,591,942 non-zero elements, data type 'float32', and first non-zero element = 1
obs: _index, bc1_well, bc1_wind, bc2_well, bc2_wind, bc3_well, bc3_wind, cell_type, cytokine, donor, gene_count, log1p_n_genes_by_counts,
     log1p_total_counts, log1p_total_counts_MT, mread_count, pct_counts_MT, sample, species, total_counts_MT, treatment, tscp_count
var: _index, n_cells
```

## Loading data

SingleCell supports reading and writing files from each of the three major single-cell ecosystems:

- scverse/Scanpy AnnData (`.h5ad`)
- Seurat (`.rds` and `.h5Seurat`)
- Bioconductor SingleCellExperiment (`.rds`)

as well as raw 10x data files (`.h5` or `.mtx`/`.mtx.gz`). See [Interoperability](interoperability.md) for details on format conversion, partial loading, and the ryp Python-R bridge.

```python
with Timer('Load dataset'):
    sc = SingleCell(
        'Parse_10M_PBMC_cytokines.h5ad', num_threads=-1,
        obs_columns=['sample', 'donor', 'cell_type', 'treatment', 'cytokine'])
```
```none
Load dataset...
Load dataset took 1m 625ms
```

`num_threads` controls parallelism for this load and all subsequent operations on the dataset. The default (`-1`) uses all available cores. On shared machines like HPC clusters, setting an explicit value avoids contention with other jobs.`num_threads` can also be set for individual steps.

A quick look at what was loaded:

```python
sc.peek_obs()
```
```none
 column     value
 _index     89_103_005__s1
 sample     Donor10_4-1BBL
 donor      Donor10
 cell_type  CD8 Naive
 treatment  cytokine
 cytokine   4-1BBL
shape: (6, 2)
```

```python
sc.peek_var()
```
```none
 column   value
 _index   TSPAN6
 n_cells  15700
shape: (2, 2)
```

## Quality control

{meth}`~single_cell.SingleCell.qc` filters low-quality cells and removes predicted doublets in a single step. The default filters are:

- **>5% mitochondrial reads**
- **<100 genes detected**
- **Zero *MALAT1* expression** -- this nuclear lncRNA is ubiquitously expressed; absence [indicates](https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2) empty droplets or cytoplasmic fragments.
- **Predicted doublets** -- scored via [cxds](https://doi.org/10.1093/bioinformatics/btz698); threshold auto-determined from simulated doublets.

If your dataset has multiple sequencing batches, specify the batch column so doublet detection runs within each batch:

```python
with Timer('Quality control'):
    sc = sc.qc(
        remove_doublets=True, batch_column='sample',
        subset=False, allow_float=True)
```
```none
Quality control...
Starting with 9,697,974 cells.
Filtering to cells with ≤5.0% mitochondrial counts...
9,443,200 cells remain after filtering to cells with ≤5.0% mitochondrial counts.
Filtering to cells with ≥100 genes detected (with non-zero count)...
9,443,200 cells remain after filtering to cells with ≥100 genes detected.
Filtering to cells with non-zero *MALAT1* expression...
9,443,163 cells remain after filtering to cells with non-zero *MALAT1* expression.
Removing predicted doublets...
8,403,723 cells remain after removing predicted doublets.
Adding a Boolean column, obs['passed_QC'], indicating which cells passed QC...
Quality control took 1m 9s
```

Each threshold is configurable:

```python
sc.qc(max_mito_fraction=0.10, min_genes=200, nonzero_MATAL1=False,
    remove_doublets=False, allow_float=True)
```

`subset=True` (default) removes failing cells but roughly doubles peak memory by copying `X`. `subset=False` keeps all cells and adds a `passed_QC` column to {attr}`~single_cell.SingleCell.obs` -- downstream methods automatically ignore flagged cells via their `QC_column` argument.

{meth}`~single_cell.SingleCell.qc` sets `uns['QCed'] = True`. Downstream methods check this flag and error if QC has not been run. For pre-cleaned data, {meth}`~single_cell.SingleCell.skip_qc` sets the flag without filtering:

```python
sc.skip_qc()
```

### Exploring QC metrics

To explore data quality before filtering -- for instance, to choose thresholds or make plots -- {meth}`~single_cell.SingleCell.qc_metrics` adds `num_counts`, `num_genes`, and `mito_fraction` columns to {attr}`~single_cell.SingleCell.obs`. This is optional; {meth}`~single_cell.SingleCell.qc` calculates its own filters internally.

```python
sc = sc.qc_metrics(allow_float=True)
sc.obs.select('num_counts', 'num_genes', 'mito_fraction').describe()
```

```none
 statistic   num_counts   num_genes    mito_fraction
 count       9.697974e6   9.697974e6   9.697974e6
 null_count  0.0          0.0          0.0
 mean        4372.856645  1941.703694  0.020779
 std         3870.176441  934.460866   0.01191
 min         436.0        399.0        0.0
 25%         2014.0       1274.0       0.012927
 50%         3320.0       1795.0       0.018277
 75%         5379.0       2417.0       0.025636
 max         70055.0      7000.0       0.149981
shape: (9, 4)
```

### Standalone doublet detection

For more control over doublet scoring without re-running QC, {meth}`~single_cell.SingleCell.find_doublets` can be called separately. It adds `doublet` and `doublet_score` columns to {attr}`~single_cell.SingleCell.obs`:

```python
sc.find_doublets(batch_column='sample')
```

## Saving

```python
sc.save('qced.h5ad', overwrite=True)
```

## Summary

| Step | Method | What it does |
|---|---|---|
| Inspect | {meth}`SingleCell.ls('file.h5ad') <single_cell.SingleCell.ls>` | Preview file structure without loading |
| Load | {meth}`SingleCell('file.h5ad') <single_cell.SingleCell.__init__>` | Read data from any supported format |
| QC metrics | {meth}`sc.qc_metrics() <single_cell.SingleCell.qc_metrics>` | Add `num_counts`, `num_genes`, `mito_fraction` to `obs` |
| Quality control | {meth}`sc.qc(remove_doublets=True) <single_cell.SingleCell.qc>` | Filter low-quality cells and doublets |
| Skip QC | {meth}`sc.skip_qc() <single_cell.SingleCell.skip_qc>` | Mark pre-cleaned data as QCed |
| Save | {meth}`sc.save('out.h5ad') <single_cell.SingleCell.save>` | Write to `.h5ad`, `.rds`, `.h5Seurat`, `.h5`, or `.mtx` |
