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


brisc bridges to R through [ryp](https://github.com/Wainberg/ryp), which uses the R arrow package to move data between Python and R, so arrow is needed for any R interop. The other R packages are feature-specific: Seurat for converting Seurat objects, limma for differential expression, and SingleCellExperiment for converting SingleCellExperiment objects.


Installing brisc with conda (`conda install -c conda-forge brisc`) sets up R, arrow, and Seurat automatically. With pip you set them up yourself: first install R — [CRAN](https://cran.r-project.org) has per-platform instructions — then add arrow and Seurat from an R session with `install.packages(c("arrow", "Seurat"))`. The optional limma and SingleCellExperiment packages are installed separately, only if you need them.


If you need differential expression, install limma:


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


If you need to convert to or from SingleCellExperiment objects, install it:


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
