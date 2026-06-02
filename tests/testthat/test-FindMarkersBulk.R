test_that("FindMarkersBulk validates input Seurat object", {
  expect_error(
    FindMarkersBulk(
      seurat = mtcars,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id"
    ),
    "'seurat' must be a Seurat object"
  )
})

test_that("FindMarkersBulk rejects missing sample metadata", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersBulk(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "missing_sample"
    ),
    "'sample_ident' not found in seurat@meta.data"
  )
})

test_that("FindMarkersBulk creates output directory before design validation", {
  seurat <- create_test_seurat()
  out_dir <- tempfile("findmarkersbulk-out-")

  expect_error(
    FindMarkersBulk(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      design_formula = ~ sample_id,
      out_dir = out_dir
    ),
    "Design formula must include 'iscluster'"
  )

  expect_true(dir.exists(out_dir))
  cleanup_test_files(out_dir)
})

test_that("FindMarkersBulk rejects invalid batch variables", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersBulk(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      batch_var = "missing_batch"
    ),
    "batch_var missing_batch not found in seurat metadata"
  )
})

test_that("FindMarkersBulk rejects invalid covariates", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersBulk(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      covariates = c("treatment", "missing_covariate")
    ),
    "Covariates not found in metadata: missing_covariate"
  )
})

test_that("FindMarkersBulk completes and writes cluster result files", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkersbulk-success-")

  result <- FindMarkersBulk(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    sample_ident = "sample_id",
    test_type = "Wald",
    expfilt_counts = 1,
    expfilt_freq = 0.25,
    alpha = 0.5,
    n_top_genes = 5,
    out_dir = out_dir
  )

  expect_type(result, "list")
  expect_true(file.exists(file.path(out_dir, "cluster_0_results.csv")))
  expect_true(file.exists(file.path(out_dir, "cluster_1_results.csv")))
  expect_true(file.exists(file.path(out_dir, "Top_markers.csv")))
  expect_true("all_results" %in% names(result))

  expect_true(all(c("0", "1") %in% names(result$all_results)))
  expect_true(nrow(result$all_results[["0"]]) > 0)
  expect_true(nrow(result$all_results[["1"]]) > 0)
  expect_true(all(c("feature", "log2FoldChange", "pct_in", "pct_out", "padj") %in% colnames(result$all_results[["0"]])))

  cluster0_csv <- read.csv(file.path(out_dir, "cluster_0_results.csv"))
  cluster1_csv <- read.csv(file.path(out_dir, "cluster_1_results.csv"))
  expect_true(nrow(cluster0_csv) > 0)
  expect_true(nrow(cluster1_csv) > 0)
  expect_true(all(c("feature", "log2FoldChange", "pct_in", "pct_out", "padj") %in% colnames(cluster0_csv)))

  top_markers <- read.csv(file.path(out_dir, "Top_markers.csv"), stringsAsFactors = FALSE)
  expect_true(nrow(top_markers) > 0)

  cleanup_test_files(out_dir)
})

test_that("FindMarkersBulk fails on non-integer pseudobulk counts", {
  seurat <- create_non_integer_test_seurat()
  out_dir <- tempfile("findmarkersbulk-noninteger-")

  expect_error(
    FindMarkersBulk(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      test_type = "Wald",
      expfilt_counts = 1,
      expfilt_freq = 0.25,
      alpha = 0.5,
      n_top_genes = 5,
      out_dir = out_dir
    ),
    "Non-integer pseudobulk counts detected"
  )

  cleanup_test_files(out_dir)
})

test_that("FindMarkersBulk LRT with batch keeps fold-change direction after shrinkage", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkersbulk-lrt-batch-")

  result <- FindMarkersBulk(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    sample_ident = "sample_id",
    batch_var = "batch",
    test_type = "LRT",
    expfilt_counts = 1,
    expfilt_freq = 0.25,
    alpha = 0.5,
    n_top_genes = 5,
    out_dir = out_dir
  )

  expect_true(all(c("0", "1") %in% names(result$all_results)))

  cluster0 <- result$all_results[["0"]]
  cluster1 <- result$all_results[["1"]]

  expect_true(all(c("log2FoldChange", "log2FoldChange_raw") %in% colnames(cluster0)))
  expect_true(all(c("log2FoldChange", "log2FoldChange_raw") %in% colnames(cluster1)))

  keep0 <- !is.na(cluster0$log2FoldChange) & !is.na(cluster0$log2FoldChange_raw)
  keep1 <- !is.na(cluster1$log2FoldChange) & !is.na(cluster1$log2FoldChange_raw)

  expect_true(any(keep0))
  expect_true(any(keep1))

  sign0 <- sign(cluster0$log2FoldChange[keep0]) == sign(cluster0$log2FoldChange_raw[keep0]) |
    cluster0$log2FoldChange[keep0] == 0 | cluster0$log2FoldChange_raw[keep0] == 0
  sign1 <- sign(cluster1$log2FoldChange[keep1]) == sign(cluster1$log2FoldChange_raw[keep1]) |
    cluster1$log2FoldChange[keep1] == 0 | cluster1$log2FoldChange_raw[keep1] == 0

  expect_true(all(sign0))
  expect_true(all(sign1))

  expect_true(any(abs(cluster0$log2FoldChange[keep0] - cluster0$log2FoldChange_raw[keep0]) > 1e-8))
  expect_true(any(abs(cluster1$log2FoldChange[keep1] - cluster1$log2FoldChange_raw[keep1]) > 1e-8))

  cleanup_test_files(out_dir)
})
