# Installation

brisc supports Linux, macOS, and Windows on Python 3.9+. Install it with conda or pip:

::::{tab-set}
:::{tab-item} conda (recommended)
```bash
conda install -c conda-forge brisc
```
:::
:::{tab-item} pip
```bash
pip install brisc
```
:::
::::

conda is recommended because it also sets up the fast MKL BLAS and the R packages brisc uses (both covered below); with pip you must handle those yourself.

## BLAS and threading

A few key steps — nearest-neighbor search, harmonization, and label transfer — rely on BLAS. conda's SciPy comes with **MKL** (fastest, no thread cap); pip's comes with **OpenBLAS**, which brisc limits to 64 threads. On a pip install with more than 64 cores, switch to MKL:

```bash
conda install "libblas=*=*mkl" scipy
```

To check which backend you have:

```python
import brisc
from threadpoolctl import threadpool_info
print(sorted({pool['internal_api'] for pool in threadpool_info()}))
```

`mkl` is the fast path; `openblas` means the 64-thread cap applies.

## R packages

brisc bridges to R through [ryp](https://github.com/Wainberg/ryp) for differential expression and conversion of Seurat, SingleCellExperiment, and `.rds` objects. conda installs the needed R packages (and R itself) automatically; with pip, you must install them yourself:

```bash
conda install -c conda-forge -c bioconda r-arrow bioconductor-limma r-seurat bioconductor-singlecellexperiment
```

Only `r-arrow` is always required, so drop any of the others you won't use.
