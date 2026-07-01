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

conda is recommended because it also sets up the fast MKL BLAS and some of the R packages brisc uses (both covered below); with pip you must handle those yourself.

## R packages

brisc's R integration is optional. You need it only for differential expression, or to work with Seurat and SingleCellExperiment data — both reading and writing their `.rds` files and converting their in-memory objects. If you do neither, you can skip this section.

The integration runs through [ryp](https://github.com/Wainberg/ryp), which moves data between Python and R via the R arrow package, so arrow is required whenever R is involved. Each feature then needs its own package: limma for differential expression, and Seurat or SingleCellExperiment for the corresponding data.

Installing brisc with conda (`conda install -c conda-forge brisc`) sets up R, arrow, and Seurat automatically. With pip you set them up yourself: first install R — [CRAN](https://cran.r-project.org) has per-platform instructions — then add arrow and Seurat from an R session with `install.packages(c("arrow", "Seurat"))`.

For differential expression, install limma:

::::{tab-set}
:::{tab-item} conda (Linux/macOS)
```bash
conda install -c bioconda bioconductor-limma
```
:::
:::{tab-item} BiocManager (all platforms)
```bash
R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("limma")'
```
:::
::::

For SingleCellExperiment data, install it:

::::{tab-set}
:::{tab-item} conda (Linux/macOS)
```bash
conda install -c bioconda bioconductor-singlecellexperiment
```
:::
:::{tab-item} BiocManager (all platforms)
```bash
R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("SingleCellExperiment")'
```
:::
::::

## BLAS and threading

A few key steps — nearest-neighbor search, harmonization, and label transfer — rely on BLAS. On machines with x86 processors (most Linux and Windows machines), conda allows SciPy to be installed with **MKL BLAS**, which is highly optimized. (You can install this manually with `conda install "libblas=*=*mkl" scipy`, although installing brisc through conda takes care of this for you.)

However, pip's SciPy comes with **OpenBLAS**, which is less optimized and only supports up to 64 threads. To check which backend you have:

```python
import brisc
from threadpoolctl import threadpool_info
print(sorted({pool['internal_api'] for pool in threadpool_info()}))
```

`mkl` is the fast path; `openblas` means the 64-thread cap applies.
