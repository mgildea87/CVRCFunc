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
