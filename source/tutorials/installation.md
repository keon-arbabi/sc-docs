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


A few key steps — nearest-neighbor search, harmonization, and label transfer — rely on BLAS. On machines with x86 processors (most Linux and Windows machines), conda allows SciPy to be installed with **MKL BLAS**, which is highly optimized. (You can install this manually with `conda install "libblas=*=*mkl" scipy`, although installing brisc through conda takes care of this for you.)


However, pip's SciPy comes with **OpenBLAS**, which is less optimized and only supports up to 64 threads. To check which backend you have:


```python
import brisc
from threadpoolctl import threadpool_info
print(sorted({pool['internal_api'] for pool in threadpool_info()}))
```


`mkl` is the fast path; `openblas` means the 64-thread cap applies.


## R packages


brisc bridges to R through [ryp](https://github.com/Wainberg/ryp) for differential expression (limma) and for converting Seurat and SingleCellExperiment objects. conda installs R and these packages automatically; with pip you set them up yourself. First install R — [CRAN](https://cran.r-project.org) has per-platform instructions — then add the R packages you need from an R session:


```r
# required by ryp for all Python-R data transfer
install.packages("arrow")


# differential expression
install.packages("BiocManager")
BiocManager::install("limma")


# convert Seurat / SingleCellExperiment objects
install.packages("Seurat")
BiocManager::install("SingleCellExperiment")
```
