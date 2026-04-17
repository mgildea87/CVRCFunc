if (!exists("BetterVlnPlot", mode = "function")) {
  pkgload::load_all(".", quiet = TRUE)
}

test_that("BetterVlnPlot returns a ggplot object", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, features = "GeneA")

  expect_s3_class(plot, "ggplot")
})

test_that("BetterVlnPlot keeps one facet row per cell-feature value", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, features = c("GeneA", "GeneB"), ncol = 2)

  expect_setequal(unique(plot$data$variable), c("GeneA", "GeneB"))
  expect_equal(nrow(plot$data), ncol(seurat) * 2)
})

test_that("BetterVlnPlot uses requested clustering column", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, clus_ident = "cluster_label", features = "GeneA")

  expect_setequal(unique(plot$data$cluster_label), c("alpha", "beta"))
})

test_that("BetterVlnPlot filters to requested identities", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, features = "GeneA", idents_to_plot = "1")

  expect_setequal(unique(plot$data$seurat_clusters), "1")
  expect_equal(nrow(plot$data), 4)
})

test_that("BetterVlnPlot includes condition metadata when splitting", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(
    seurat,
    features = c("GeneA", "GeneB"),
    condition_ident = "treatment"
  )

  expect_true("treatment" %in% colnames(plot$data))
  expect_setequal(unique(plot$data$treatment), c("ctrl", "stim"))
})

test_that("BetterVlnPlot omits quasirandom layer when dot_size is zero", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, features = "GeneA", dot_size = 0)

  expect_length(plot$layers, 1)
  expect_match(class(plot$layers[[1]]$geom)[1], "GeomViolin")
})

test_that("BetterVlnPlot applies the requested y axis title", {
  seurat <- create_test_seurat()

  plot <- BetterVlnPlot(seurat, features = "GeneA", y_axis_title = "custom title")

  expect_equal(plot$labels$y, "custom title")
})

test_that("BetterVlnPlot errors clearly for missing features", {
  seurat <- create_test_seurat()

  expect_error(
    BetterVlnPlot(seurat, features = "MissingGene"),
    "None of the requested features were found in the assay data"
  )
})
