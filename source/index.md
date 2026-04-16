<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="_static/images/runner_title_wide_light.svg" alt="brisc" class="only-light" style="max-width: 80%;">
  <img src="_static/images/runner_title_wide_dark.svg" alt="brisc" class="only-dark" style="max-width: 80%;">
</div>

brisc is a high-performance library for analyzing single-cell data at scale. It prioritizes running as fast as possible on multithreaded CPU systems, strict reproducibility and a clean, user-friendly interface. On datasets of 1 to 20 million cells, it cuts the runtime of common workflows by over 60x, from hours to minutes.

## Key features

-   **Blazing fast**: Achieved through ground-up optimization of core algorithms and effective parallelism.

-   **Deterministic parallelism:** Every step produces identical floating-point results, regardless of the number of threads used.

-   **User friendly**: Sensible defaults, strict type-checking, and specific error messages keep the interface approachable.

-   **Seamless interoperability**: Natively read and write `.h5ad`, `.rds`, `.h5Seurat`, and 10x files. Convert between SingleCell, AnnData, Seurat, and SingleCellExperiment objects in memory via the [ryp](https://github.com/Wainberg/ryp) Python-R bridge.

-   **Complete toolkit**: Covers every major single-cell task — preprocessing, dimensionality reduction, harmony, label transfer, clustering, embedding, and pseudobulk differential expression.

-   **Memory efficient**: Tracks QC flagged cells instead of subsetting, avoiding a 2x peak in memory.

<div class="showcase-row">
<div>
<div class="showcase-label">Basic Workflow</div>

```python
from brisc import SingleCell

sc = SingleCell('data.h5ad')\
  .qc()\
  .hvg(batch_column='sample')\
  .normalize()\
  .pca()\
  .neighbors()\
  .shared_neighbors()\
  .cluster(resolution=[0.5, 1.0, 2.0])\
  .pacmap()
```

</div>
<div>
<div class="showcase-label">Label Transfer</div>

```python
from brisc import SingleCell

sc_ref = SingleCell('data_ref.h5ad').qc()
sc_query = SingleCell('data_query.h5ad').qc()
sc_ref, sc_query = sc_ref.hvg(sc_query)
sc_ref = sc_ref.normalize()
sc_query = sc_query.normalize()
sc_ref, sc_query = sc_ref.pca(sc_query)
sc_ref, sc_query = sc_ref.harmonize(sc_query)
sc_query = sc_query.label_transfer_from(
  sc_ref, 'cell_type',
  cell_type_column='cell_type_transferred')
```

</div>
<div>
<div class="showcase-label">Pseudobulk DE</div>

```python
from brisc import SingleCell, Pseudobulk, DE

sc = SingleCell('data.h5ad').qc()
pb = sc.pseudobulk('sample', 'cell_type')
pb = pb.qc('condition')
formula = '~ condition + sex + age_at_death'
de = pb.DE(
  formula, group='condition',
  categorical_columns=['condition', 'sex'])
```

</div>
</div>

<div class="showcase-chart">
<div class="showcase-chart-header">
<div class="showcase-chart-title">Performance · 10M cells · Parse Biosciences PBMC</div>
<div class="showcase-chart-subtitle">192 CPUs, 755 GB RAM</div>
</div>
<div id="benchmark-chart"></div>
</div>

<script>
document.body.classList.add('landing-page');
</script>

:::{toctree}
:caption: Tutorials
:hidden:
:maxdepth: 2

tutorials/index

:::

:::{toctree}
:caption: API
:hidden:
:maxdepth: 4

api/single_cell/index
api/pseudobulk/index
api/de/index

:::
