<div style="text-align: center; margin-bottom: 1.5rem;">
  <img src="_static/images/runner_title_light.svg" alt="brisc" class="only-light" style="max-width: 50%;">
  <img src="_static/images/runner_title_dark.svg" alt="brisc" class="only-dark" style="max-width: 50%;">
</div>

brisc is a high-performance library for analyzing single-cell data at scale. It prioritizes running as fast as possible on multithreaded CPU systems, while prioritizing strict reproducibilty and a clean, user-friendly interface. On datasets of 1 to 20 million cells, it cuts the runtime of common workflows by over 60x, from hours to minutes.

Built on [Polars](https://pola.rs/) and cython, it provides a blazing fast, memory-efficient pipeline that interoperates seamlessly with the Scverse and Seurat ecosystems.

## Key features

-   **Blazing fast**: This is achieved through a combination of extensive optimization of core algorithms and  parallelism.

A multi-threaded reader loads massive HDf5 in seconds to minutes.

-   **Deterministic Parallelism** Every step is specifically designed to give identical float-point results, regardless of the number of threads used.

-   **Seamless Interoperability**: Natively read and write Seurat (`.rds`, `.h5Seurat`) and AnnData (`.h5ad`) files, bridging the Python and R ecosystems without intermediate conversion steps.

-   **Complete Pipeline**: A unified toolkit spanning QC (including doublet detection), batch correction, clustering, and pseudobulk differential expression via limma-voom.

-   **Memory Efficient**: brisc tracks QCed cells rather than subsetting to reduce peak memory useage,

-   **User Friendly**: Designed to be beginner-friendly. Detailed, specific error messsages. Type-checking. Sensible defaults. 


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
