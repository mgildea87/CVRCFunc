test_that("FindMarkers rejects invalid clus_ident", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkers(
      seurat = seurat,
      clus_ident = "missing_column",
      group_1 = "0",
      group_2 = "1",
      sample_ident = "sample_id"
    ),
    "clus_ident missing_column not found in metadata"
  )
})

test_that("FindMarkers rejects invalid group identities", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkers(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      group_1 = "bad_group",
      group_2 = "1",
      sample_ident = "sample_id"
    ),
    "One or both groups not found in the specified identity column"
  )
})

test_that("FindMarkers creates output directory before later validation errors", {
  seurat <- create_test_seurat()
  out_dir <- tempfile("findmarkers-out-")

  expect_error(
    FindMarkers(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      group_1 = "0",
      group_2 = "1",
      sample_ident = "sample_id",
      batch_var = "missing_batch",
      out_dir = out_dir
    ),
    "batch_var missing_batch not found in seurat metadata"
  )

  expect_true(dir.exists(out_dir))
  cleanup_test_files(out_dir)
})

test_that("FindMarkers rejects missing batch variables", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkers(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      group_1 = "0",
      group_2 = "1",
      sample_ident = "sample_id",
      batch_var = "missing_batch"
    ),
    "batch_var missing_batch not found in seurat metadata"
  )
})

test_that("FindMarkers rejects missing covariates", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkers(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      group_1 = "0",
      group_2 = "1",
      sample_ident = "sample_id",
      covariates = c("treatment", "missing_covariate")
    ),
    "Covariates not found in metadata: missing_covariate"
  )
})

test_that("FindMarkers completes and writes expected result files", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkers-success-")

  result <- FindMarkers(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    group_1 = "0",
    group_2 = "1",
    sample_ident = "sample_id",
    test_type = "Wald",
    expfilt_counts = 1,
    expfilt_freq = 0.25,
    alpha = 0.5,
    out_dir = out_dir
  )

  expect_type(result, "list")
  expect_true(file.exists(file.path(out_dir, "results.csv")))
  expect_true(file.exists(file.path(out_dir, "diagnostic_plots.pdf")))
  expect_true("results" %in% names(result))
  expect_true("dds" %in% names(result))

  expect_true(nrow(result$results) > 0)
  expect_true(all(c("feature", "log2FoldChange", "pct_in", "pct_out", "padj") %in% colnames(result$results)))

  csv_results <- read.csv(file.path(out_dir, "results.csv"))
  expect_true(nrow(csv_results) > 0)
  expect_true(all(c("feature", "log2FoldChange", "pct_in", "pct_out", "padj") %in% colnames(csv_results)))

  cluster0_genes <- paste0("Gene", sprintf("%03d", 1:40))
  cluster1_genes <- paste0("Gene", sprintf("%03d", 41:80))
  mean_lfc_cluster0 <- mean(result$results$log2FoldChange[result$results$feature %in% cluster0_genes], na.rm = TRUE)
  mean_lfc_cluster1 <- mean(result$results$log2FoldChange[result$results$feature %in% cluster1_genes], na.rm = TRUE)
  expect_lt(mean_lfc_cluster0, 0)
  expect_gt(mean_lfc_cluster1, 0)

  cleanup_test_files(out_dir)
})
