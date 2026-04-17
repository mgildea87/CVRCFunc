# Test setup and fixtures for CVRCFunc

create_test_seurat <- function() {
  counts <- Matrix::Matrix(matrix(
    c(
      12, 11, 10,  9,  1,  1,  0,  0,
       2,  3,  2,  1, 11, 12, 10,  9,
       6,  5,  7,  6,  5,  4,  6,  5,
       1,  0,  1,  0,  8,  7,  9,  8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(
      c("GeneA", "GeneB", "GeneC", "GeneD"),
      paste0("Cell", seq_len(8))
    )
  ), sparse = TRUE)

  meta.data <- data.frame(
    seurat_clusters = factor(c("0", "0", "0", "0", "1", "1", "1", "1")),
    cluster_label = factor(c("alpha", "alpha", "alpha", "alpha", "beta", "beta", "beta", "beta")),
    treatment = factor(c("ctrl", "stim", "ctrl", "stim", "ctrl", "stim", "ctrl", "stim")),
    sample_id = factor(c("S1", "S2", "S3", "S4", "S1", "S2", "S3", "S4")),
    row.names = colnames(counts)
  )

  seurat <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta.data)
  seurat <- Seurat::NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seurat
}

create_de_test_seurat <- function() {
  set.seed(42)

  genes <- paste0("Gene", sprintf("%03d", seq_len(140)))
  samples <- c("ctrl_1", "ctrl_2", "stim_1", "stim_2")
  clusters <- c("0", "1")
  cells_per_sample_cluster <- 2

  meta_rows <- expand.grid(
    sample_id = samples,
    seurat_clusters = clusters,
    rep = seq_len(cells_per_sample_cluster),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  meta_rows$cell <- paste0("Cell", seq_len(nrow(meta_rows)))
  meta_rows$treatment <- ifelse(grepl("^ctrl", meta_rows$sample_id), "ctrl", "stim")
  meta_rows$batch <- ifelse(meta_rows$sample_id %in% c("ctrl_1", "stim_1"), "B1", "B2")
  meta_rows$cluster_label <- ifelse(meta_rows$seurat_clusters == "0", "alpha", "beta")

  counts <- matrix(
    rpois(length(genes) * nrow(meta_rows), lambda = 8),
    nrow = length(genes),
    dimnames = list(genes, meta_rows$cell)
  )

  cluster0_cells <- meta_rows$cell[meta_rows$seurat_clusters == "0"]
  cluster1_cells <- meta_rows$cell[meta_rows$seurat_clusters == "1"]
  stim_cells <- meta_rows$cell[meta_rows$treatment == "stim"]
  ctrl_cells <- meta_rows$cell[meta_rows$treatment == "ctrl"]

  counts[1:40, cluster0_cells] <- counts[1:40, cluster0_cells] + 60L
  counts[41:80, cluster1_cells] <- counts[41:80, cluster1_cells] + 60L
  counts[81:110, intersect(cluster0_cells, stim_cells)] <- counts[81:110, intersect(cluster0_cells, stim_cells)] + 50L
  counts[111:140, intersect(cluster1_cells, ctrl_cells)] <- counts[111:140, intersect(cluster1_cells, ctrl_cells)] + 50L

  meta.data <- data.frame(
    seurat_clusters = factor(meta_rows$seurat_clusters),
    cluster_label = factor(meta_rows$cluster_label),
    treatment = factor(meta_rows$treatment),
    sample_id = factor(meta_rows$sample_id),
    batch = factor(meta_rows$batch),
    row.names = meta_rows$cell
  )

  seurat <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    meta.data = meta.data
  )
  seurat <- Seurat::NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  seurat
}

cleanup_test_files <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE)
  }
}
