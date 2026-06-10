# Tutorials

These tutorials all use the same Parse Biosciences ~10 million cell PBMC dataset. New to brisc? Start with **Installation**, then **Basic Workflow**.

| Tutorial | Description |
|---|---|
| [Installation](installation.md) | Install brisc and its optional R dependencies. |
| [Basic Workflow](basic_workflow.md) | A standard analysis from raw counts to clusters and markers. |
| [Integration and Label Transfer](integration_and_label_transfer.md) | Map an annotated reference onto a query and transfer its cell-type labels. |
| [Differential Expression](differential_expression.md) | Pseudobulk differential expression between conditions with limma-voom. |
| [Interoperability](interoperability.md) | Move data between brisc, scanpy/AnnData, Seurat, and SingleCellExperiment. |
| [Data Manipulation](data_manipulation.md) | A polars cheatsheet for inspecting, subsetting, and reshaping a dataset. |

:::{toctree}
:hidden:

installation
basic_workflow
integration_and_label_transfer
differential_expression
interoperability
data_manipulation
:::
