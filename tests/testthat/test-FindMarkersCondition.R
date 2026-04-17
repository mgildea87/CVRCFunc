test_that("FindMarkersCondition validates condition identity parameter", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "missing_condition",
      conditions = c("ctrl", "stim")
    ),
    "'condition_ident' missing_condition not found in metadata"
  )
})

test_that("FindMarkersCondition requires two explicit conditions for Wald tests", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "treatment",
      conditions = NULL,
      test_type = "Wald"
    ),
    "conditions must be specified when test_type = 'Wald'"
  )
})

test_that("FindMarkersCondition rejects unknown conditions", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "treatment",
      conditions = c("ctrl", "missing_condition")
    ),
    "One or both conditions not found in the data"
  )
})

test_that("FindMarkersCondition creates output directory before design validation", {
  seurat <- create_test_seurat()
  out_dir <- tempfile("findmarkerscondition-out-")

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "treatment",
      conditions = c("ctrl", "stim"),
      design_formula = ~ sample_id,
      out_dir = out_dir
    ),
    "Design formula must include 'treatment'"
  )

  expect_true(dir.exists(out_dir))
  cleanup_test_files(out_dir)
})

test_that("FindMarkersCondition rejects invalid batch variables", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "treatment",
      conditions = c("ctrl", "stim"),
      batch_var = "missing_batch"
    ),
    "batch_var missing_batch not found in seurat metadata"
  )
})

test_that("FindMarkersCondition rejects invalid covariates", {
  seurat <- create_test_seurat()

  expect_error(
    FindMarkersCondition(
      seurat = seurat,
      clus_ident = "seurat_clusters",
      sample_ident = "sample_id",
      condition_ident = "treatment",
      conditions = c("ctrl", "stim"),
      covariates = c("cluster_label", "missing_covariate")
    ),
    "Covariates not found in metadata: missing_covariate"
  )
})

test_that("FindMarkersCondition completes and writes summary outputs", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkerscondition-success-")

  result <- FindMarkersCondition(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    sample_ident = "sample_id",
    condition_ident = "treatment",
    conditions = c("stim", "ctrl"),
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
  expect_true(file.exists(file.path(out_dir, "Summary_results.csv")))
  expect_true("summary" %in% names(result))

  expect_true(nrow(result$summary) == 2)
  expect_true(all(c("cluster", "sig_up_in_stim", "sig_down_in_stim") %in% colnames(result$summary)))

  summary_tbl <- read.csv(file.path(out_dir, "Summary_results.csv"))
  expect_true(nrow(summary_tbl) == 2)
  expect_true(all(c("cluster", "sig_up_in_stim", "sig_down_in_stim") %in% colnames(summary_tbl)))

  row0 <- summary_tbl[summary_tbl$cluster == 0 | summary_tbl$cluster == "0", , drop = FALSE]
  row1 <- summary_tbl[summary_tbl$cluster == 1 | summary_tbl$cluster == "1", , drop = FALSE]
  expect_gt(row0$sig_up_in_stim, row1$sig_up_in_stim)
  expect_gt(row1$sig_down_in_stim, row0$sig_down_in_stim)

  cluster0_csv <- read.csv(file.path(out_dir, "cluster_0_results.csv"))
  expect_true(all(c("feature", "log2FoldChange", "pct_in_stim", "pct_in_ctrl", "padj") %in% colnames(cluster0_csv)))

  cleanup_test_files(out_dir)
})
