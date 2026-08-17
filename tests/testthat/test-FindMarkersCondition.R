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

test_that("FindMarkersCondition removes genes below pct.in before DESeq", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkerscondition-pctin-")

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
    pct.in = 0.5,
    out_dir = out_dir
  )

  expect_type(result, "list")
  expect_true(all(vapply(result$all_results, function(x) {
    all(x$pct_in_stim >= 0.5 | is.na(x$pct_in_stim))
  }, logical(1))))
  expect_true(all(vapply(result$all_results, function(x) {
    all(x$pct_out_stim >= 0.5 | is.na(x$pct_out_stim))
  }, logical(1))))

  cleanup_test_files(out_dir)
})

test_that("FindMarkersCondition writes console output to a log file by default", {
  seurat <- create_de_test_seurat()
  out_dir <- tempfile("findmarkerscondition-log-")

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
  expect_true(file.exists(file.path(out_dir, "analysis.log")))
  log_contents <- readLines(file.path(out_dir, "analysis.log"), warn = FALSE)
  expect_true(any(grepl("Processing cluster", log_contents, fixed = TRUE)))

  cleanup_test_files(out_dir)
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

test_that("FindMarkersCondition fails on non-integer pseudobulk counts", {
  seurat <- create_non_integer_test_seurat()
  out_dir <- tempfile("findmarkerscondition-noninteger-")

  expect_error(
    FindMarkersCondition(
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
    ),
    "Non-integer pseudobulk counts detected"
  )

  cleanup_test_files(out_dir)
})

test_that("FindMarkersCondition LRT with batch uses supplied condition direction", {
  seurat <- create_de_test_seurat()
  seurat$treatment <- factor(seurat$treatment, levels = c("stim", "ctrl"))
  out_dir <- tempfile("findmarkerscondition-lrt-batch-")

  result <- FindMarkersCondition(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    sample_ident = "sample_id",
    condition_ident = "treatment",
    conditions = c("stim", "ctrl"),
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

    expect_true("log2FoldChange" %in% colnames(cluster0))
    expect_true("log2FoldChange" %in% colnames(cluster1))
    expect_false("log2FoldChange_raw" %in% colnames(cluster0))
    expect_false("log2FoldChange_raw" %in% colnames(cluster1))
    expect_true("log2FoldChange_treatment_stim_vs_ctrl" %in% colnames(cluster0))
    expect_true("log2FoldChange_treatment_stim_vs_ctrl" %in% colnames(cluster1))
    expect_equal(cluster0$log2FoldChange, cluster0$log2FoldChange_treatment_stim_vs_ctrl)
    expect_equal(cluster1$log2FoldChange, cluster1$log2FoldChange_treatment_stim_vs_ctrl)

  genes_up_in_stim_cluster0 <- paste0("Gene", sprintf("%03d", 81:110))
  genes_up_in_ctrl_cluster1 <- paste0("Gene", sprintf("%03d", 111:140))

  mean_lfc_cluster0 <- mean(
    cluster0$log2FoldChange[cluster0$feature %in% genes_up_in_stim_cluster0],
    na.rm = TRUE
  )
  mean_lfc_cluster1 <- mean(
    cluster1$log2FoldChange[cluster1$feature %in% genes_up_in_ctrl_cluster1],
    na.rm = TRUE
  )

  expect_gt(mean_lfc_cluster0, 0)
  expect_lt(mean_lfc_cluster1, 0)

  cleanup_test_files(out_dir)
})

test_that("FindMarkersCondition LRT with >3 condition levels respects supplied contrast", {
  set.seed(42)

  genes <- paste0("Gene", sprintf("%03d", seq_len(140)))
  samples <- c(
    "ctrl_1", "ctrl_2",
    "other1_1", "other1_2",
    "stim_1", "stim_2",
    "other2_1", "other2_2"
  )
  condition_levels <- c("ctrl", "other1", "stim", "other2")
  sample_conditions <- c("ctrl", "ctrl", "other1", "other1", "stim", "stim", "other2", "other2")

  meta_rows <- expand.grid(
    sample_id = samples,
    seurat_clusters = c("0", "1"),
    rep = seq_len(2),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  meta_rows$cell <- paste0("Cell", seq_len(nrow(meta_rows)))
  meta_rows$treatment <- sample_conditions[match(meta_rows$sample_id, samples)]

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
    cluster_label = factor(ifelse(meta_rows$seurat_clusters == "0", "alpha", "beta")),
    treatment = factor(meta_rows$treatment, levels = condition_levels),
    sample_id = factor(meta_rows$sample_id, levels = samples),
    row.names = meta_rows$cell
  )

  seurat <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    meta.data = meta.data
  )
  seurat <- Seurat::NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

  out_dir <- tempfile("findmarkerscondition-lrt-multilevel-")
  result <- FindMarkersCondition(
    seurat = seurat,
    clus_ident = "seurat_clusters",
    sample_ident = "sample_id",
    condition_ident = "treatment",
    conditions = c("stim", "ctrl"),
    test_type = "LRT",
    expfilt_counts = 1,
    expfilt_freq = 0.25,
    alpha = 0.5,
    n_top_genes = 5,
    out_dir = out_dir
  )

  cluster0 <- result$all_results[["0"]]
  cluster1 <- result$all_results[["1"]]

  lfc_cols_cluster0 <- grep("^log2FoldChange_treatment_", colnames(cluster0), value = TRUE)
  lfc_cols_cluster1 <- grep("^log2FoldChange_treatment_", colnames(cluster1), value = TRUE)

  mean_lfc_cluster0 <- mean(
    cluster0$log2FoldChange[cluster0$feature %in% paste0("Gene", sprintf("%03d", 81:110))],
    na.rm = TRUE
  )
  mean_lfc_cluster1 <- mean(
    cluster1$log2FoldChange[cluster1$feature %in% paste0("Gene", sprintf("%03d", 111:140))],
    na.rm = TRUE
  )

  expect_gt(mean_lfc_cluster0, 0)
  expect_lt(mean_lfc_cluster1, 0)
  expect_true("log2FoldChange_treatment_stim_vs_ctrl" %in% lfc_cols_cluster0)
  expect_true("log2FoldChange_treatment_stim_vs_ctrl" %in% lfc_cols_cluster1)
  expect_equal(length(lfc_cols_cluster0), 3L)
  expect_equal(length(lfc_cols_cluster1), 3L)
  expect_equal(cluster0$log2FoldChange, cluster0$log2FoldChange_treatment_stim_vs_ctrl)
  expect_equal(cluster1$log2FoldChange, cluster1$log2FoldChange_treatment_stim_vs_ctrl)

  cleanup_test_files(out_dir)
})
