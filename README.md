# CVRCFunc

CVRCFunc is a lightweight R package for CVRC single-cell analysis helpers. It contains pseudobulk differential expression wrappers, clustering sweep utilities, and plotting helpers built around `Seurat`, `DESeq2`, and `ggplot2`.

## Installation

Install from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("mgildea87/CVRCFunc")
```

Install locally from source:

```r
devtools::install_local("/path/to/CVRCFunc")
```

## Main functions

### `FindMarkersBulk()`
Performs pseudobulk differential expression for every cluster using sample-level aggregation and `DESeq2`.
- `seurat`: a `Seurat` object.
- `clus_ident`: cluster identity column in `seurat@meta.data`.
- `sample_ident`: sample identifier column used for pseudobulking.
- optional `batch_var`, `covariates`, or custom `design_formula`.
- outputs CSV and PDF QC files into `out_dir`.

### `FindMarkersCondition()`
Runs pseudobulk DE per cluster between two conditions.
- `condition_ident`: metadata column defining condition labels.
- `conditions`: a length-two vector of condition values, e.g. `c("treated", "control")`.
- positive `log2FC` means higher expression in `conditions[1]`.
- supports `batch_var`, `covariates`, and custom `design_formula`.

### `FindMarkers()`
Compares two groups of cells defined by a metadata identity.
- `group_1` and `group_2` are values from `clus_ident`.
- `sample_ident` controls pseudobulk aggregation.
- useful for pairwise cluster/group comparisons.

### `FindClusterSweep()`
Runs `Seurat::FindClusters()` across a vector of resolutions and generates cluster QC plots.
- requires a neighbor graph already computed (run `FindNeighbors()` first).
- returns the input `Seurat` object with new clustering columns.
- generates a PDF with cluster counts, silhouette plots, UMAP coloring, and modularity diagnostics.

### `BetterVlnPlot()`
Creates refined violin plots for feature expression from a `Seurat` object.
- supports optional `condition_ident` to split by condition.
- accepts `idents_to_plot`, `ncol`, `dot_size`, and `y_axis_title`.
- returns a `ggplot2` object.

## Example usage

```r
library(CVRCFunc)

# Bulk cluster DE
bulk_res <- FindMarkersBulk(
  seurat = seu,
  clus_ident = "seurat_clusters",
  sample_ident = "sample_id",
  batch_var = "batch",
  covariates = c("age", "sex"),
  out_dir = "FindMarkersBulk_output"
)

# Condition DE within clusters
cond_res <- FindMarkersCondition(
  seurat = seu,
  clus_ident = "seurat_clusters",
  sample_ident = "sample_id",
  condition_ident = "condition",
  conditions = c("treated", "control"),
  out_dir = "FindMarkersCondition_output"
)

# Pairwise group comparison
pair_res <- FindMarkers(
  seurat = seu,
  clus_ident = "seurat_clusters",
  group_1 = "0",
  group_2 = "1",
  sample_ident = "sample_id",
  batch_var = "batch",
  out_dir = "FindMarkers_output"
)

# Cluster resolution sweep
seu <- FindClusterSweep(
  seurat = seu,
  assay = "RNA",
  resolutions = c(0.2, 0.4, 0.6, 0.8),
  algorithm = 1,
  reduction = "pca",
  plot_reduction = "umap",
  file_name = "cluster_sweep"
)

# Better violin plot
p <- BetterVlnPlot(
  seurat = seu,
  clus_ident = "seurat_clusters",
  features = c("CD3D", "GNLY"),
  assay = "RNA",
  condition_ident = "treatment",
  ncol = 2,
  dot_size = 0.5
)
print(p)
```

## Notes

- Requires `Seurat`, `DESeq2`, `ggplot2`, `pheatmap`, `tidyr`, `ggbeeswarm`, and related packages.
- `FindClusterSweep()` assumes a neighbor graph already exists in the `Seurat` object.
- Functions create output directories automatically when needed.
