#' Differential expression testing between 2 specified groups of cells via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param group_1 Identity of cells within \code{clus_ident} to compare
#' @param group_2 Identity of cells within \code{clus_ident} to compare
#' @param expfilt_freq genes that have greater than \code{expfilt_counts} in greater than \code{expfilt_freq} fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @param expfilt_counts genes with less than \code{expfilt_counts} in \code{expfilt_freq * sample number} will be removed from DESeq2 model. 1 by default.
#' @param alpha FDR adjusted p-value threshold for significance in plotting. 0.1 by default.
#' @param out_dir Name of output directory
#' @param assay Which assay to use. RNA by default. I added this parameter to enable use of ADT data when desired.
#' @return .csv files with marker genes per \code{clus_ident}. .pdf files with diagnostic plots
#' @import Seurat pheatmap DESeq2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom magrittr set_colnames
#' @importFrom tidyr pivot_longer
#' @export

FindMarkers <- function(seurat,
                        clus_ident,
                        group_1,
                        group_2,
                        sample_ident,
                        batch_var = NULL,
                        covariates = NULL,
                        design_formula = NULL,
                        expfilt_counts = 1,
                        expfilt_freq = 0.5,
                        out_dir = "FindMarkers",
                        alpha = 0.1,
                        assay = 'RNA',
                        test_type = "LRT",
                        contrast_level = NULL){

  start <- Sys.time()
  coef <- variable <- value <- NULL

  # Create output directory
  dir.create(out_dir, showWarnings = FALSE)

  # Validate inputs
  if (!clus_ident %in% colnames(seurat@meta.data)) {
    stop(paste("clus_ident", clus_ident, "not found in metadata"))
  }

  if (!all(c(group_1, group_2) %in% unique(seurat@meta.data[[clus_ident]]))) {
    stop("One or both groups not found in the specified identity column")
  }

  if (!is.null(batch_var) && !batch_var %in% colnames(seurat@meta.data)) {
    stop(paste("batch_var", batch_var, "not found in seurat metadata"))
  }

  if (!is.null(covariates)) {
    missing_covs <- covariates[!covariates %in% colnames(seurat@meta.data)]
    if (length(missing_covs) > 0) {
      stop(paste("Covariates not found in metadata:", paste(missing_covs, collapse = ", ")))
    }
  }

  # Set identities and subset
  Idents(seurat) <- clus_ident
  seurat <- subset(seurat, idents = c(group_1, group_2))
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  cat("=== FindMarkers Analysis ===\n")
  cat(paste("Comparing:", group_1, "vs", group_2, "\n"))
  cat(paste("Start time:", Sys.time(), "\n\n"))

  # Drop unused factor levels
  seurat@meta.data[, clus_ident] <- droplevels(factor(seurat@meta.data[, clus_ident]))

  # Build design formula
  if (is.null(design_formula)) {
    design_terms <- "cluster"

    if (!is.null(batch_var)) {
      design_terms <- c(design_terms, batch_var)
      cat(paste("Adding batch variable to design:", batch_var, "\n"))
    }

    if (!is.null(covariates)) {
      design_terms <- c(design_terms, covariates)
      cat(paste("Adding covariates to design:", paste(covariates, collapse = ", "), "\n"))
    }

    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))

    # Create reduced formula for LRT
    if (test_type == "LRT") {
      reduced_terms <- design_terms[design_terms != "cluster"]
      if (length(reduced_terms) > 0) {
        reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
      } else {
        reduced_formula <- as.formula("~ 1")
      }
    }
  } else {
    # User provided custom formula
    cat("Using custom design formula\n")
    if (test_type == "LRT") {
      # Try to create reduced formula by removing cluster
      reduced_formula <- update(design_formula, ~ . - cluster)
      # If reduced formula is empty, use intercept only
      if (length(attr(terms(reduced_formula), "term.labels")) == 0) {
        reduced_formula <- as.formula("~ 1")
      }
    }
  }

  cat(paste("\nFull design formula:", deparse(design_formula), "\n"))
  if (test_type == "LRT") {
    cat(paste("Reduced design formula:", deparse(reduced_formula), "\n"))
  }
  cat(paste("Test type:", test_type, "\n\n"))

  # Print out the table of cells in each cluster-sample group
  pdf(paste(out_dir, '/', group_1, "_", group_2, "_", 'cells_per_group_HM.pdf', sep = ''))
  pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, sample_ident]),
           display_numbers = TRUE,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize_number = 8,
           main = paste("Cells per group and sample\n", group_1, "vs", group_2))
  dev.off()

  # If batch variable provided, show batch distribution
  if (!is.null(batch_var)) {
    pdf(paste(out_dir, '/', group_1, "_", group_2, "_", 'cells_per_batch_HM.pdf', sep = ''))
    pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, batch_var]),
             display_numbers = TRUE,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             fontsize_number = 8,
             main = paste("Cells per group and batch\n", group_1, "vs", group_2))
    dev.off()
  }

  # Create explicit cluster_sample identifier
  seurat$cluster_sample <- paste(seurat@meta.data[[clus_ident]],
                                 seurat@meta.data[[sample_ident]],
                                 sep = "_")

  # Pseudobulk aggregation
  cat("Aggregating to pseudobulk...\n")
  if (length(grep(seurat@version, pattern = '^4.')) == 1) {
    pb <- aggregate.Matrix(t(seurat@assays$RNA@counts),
                           groupings = seurat$cluster_sample,
                           fun = "sum")
  } else if (length(grep(seurat@version, pattern = '^5.')) == 1) {
    pb <- aggregate.Matrix(t(seurat@assays[[assay]]@layers$counts),
                           groupings = seurat$cluster_sample,
                           fun = "sum")
    colnames(pb) <- row.names(seurat@assays[[assay]])
  }

  # Create cluster counts matrix
  cluster_counts <- as.data.frame(t(as.matrix(pb)))

  # Build metadata
  cluster_metadata <- data.frame(
    cluster = sub("_.*", "", colnames(cluster_counts)),
    sample = sub("^[^_]*_", "", colnames(cluster_counts)),
    stringsAsFactors = FALSE
  )

  # Add batch and covariate information
  if (!is.null(batch_var) || !is.null(covariates)) {
    variables_to_add <- c(batch_var, covariates)
    sample_info <- unique(seurat@meta.data[, c(sample_ident, variables_to_add), drop = FALSE])
    cluster_metadata <- merge(cluster_metadata, sample_info,
                              by.x = "sample", by.y = sample_ident,
                              sort = FALSE)

    # Reorder to match original order
    cluster_metadata <- cluster_metadata[match(colnames(cluster_counts),
                                               paste(cluster_metadata$cluster,
                                                     cluster_metadata$sample, sep = "_")), ]
  }

  # Set cluster factor levels (reference is group_1)
  cluster_metadata$cluster <- factor(cluster_metadata$cluster,
                                     levels = c(group_1, group_2))

  # Convert batch and covariates to factors
  if (!is.null(batch_var)) {
    cluster_metadata[[batch_var]] <- as.factor(cluster_metadata[[batch_var]])
  }
  if (!is.null(covariates)) {
    for (cov in covariates) {
      if (is.character(cluster_metadata[[cov]])) {
        cluster_metadata[[cov]] <- as.factor(cluster_metadata[[cov]])
      }
    }
  }

  # Check sample sizes
  n_group1 <- sum(cluster_metadata$cluster == group_1)
  n_group2 <- sum(cluster_metadata$cluster == group_2)

  cat(paste("Samples -", group_1, ":", n_group1, "|", group_2, ":", n_group2, "\n"))

  if (n_group1 < 2 || n_group2 < 2) {
    stop("Each group must have at least 2 samples")
  }

  # Check for batch confounding if applicable
  if (!is.null(batch_var)) {
    batch_table <- table(cluster_metadata$cluster, cluster_metadata[[batch_var]])
    cat("\nBatch distribution:\n")
    print(batch_table)

    # Chi-square test for independence
    if (all(dim(batch_table) >= 2) && all(batch_table >= 0)) {
      chi_test <- tryCatch(chisq.test(batch_table), error = function(e) NULL)
      if (!is.null(chi_test)) {
        cat(paste("Chi-square test for batch-group independence: p =",
                  format.pval(chi_test$p.value), "\n"))
        if (chi_test$p.value < 0.05) {
          warning("Batch may be confounded with group comparison!")
        }
      }
    }

    # Warn if batches are perfectly confounded
    if (any(rowSums(batch_table > 0) == 1)) {
      warning("Groups may be confounded with batch variable!")
    }
    cat("\n")
  }

  # Create DESeq2 object
  cat("Creating DESeq2 object...\n")
  dds <- DESeqDataSetFromMatrix(cluster_counts,
                                colData = cluster_metadata,
                                design = design_formula)

  # Filter data
  cat(paste("Filtering genes: min", expfilt_counts, "counts in",
            expfilt_freq * 100, "% of samples\n"))
  keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq * nrow(cluster_metadata)
  cat(paste("Genes before filtering:", nrow(dds), "\n"))
  dds <- dds[keep, ]
  cat(paste("Genes after filtering:", nrow(dds), "\n\n"))

  if (nrow(dds) < 100) {
    stop("Too few genes remaining after filtering")
  }

  # VST transformation
  cat("Performing variance stabilizing transformation...\n")
  vst <- varianceStabilizingTransformation(dds, blind = FALSE)

  # DESeq2
  cat(paste("Running DESeq2 with", test_type, "test...\n"))
  if (test_type == "LRT") {
    dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
    res <- results(dds, alpha = alpha)
  } else if (test_type == "Wald") {
    dds <- DESeq(dds, test = "Wald")
    # Extract results for cluster coefficient
    contrast_names <- grep("cluster", resultsNames(dds), value = TRUE)
    if (length(contrast_names) == 0) {
      stop("No 'cluster' coefficient found in resultsNames(dds)")
    }
    contrast_name <- contrast_names[1]
    if (length(contrast_names) > 1) {
      message("Multiple 'cluster' coefficients found; using: ", contrast_name)
    }
    res <- results(dds, name = contrast_name, alpha = alpha)
  }

  cat(paste("Significant genes (padj <", alpha, "):", sum(res$padj < alpha, na.rm = TRUE), "\n\n"))

  # Shrink the log2 fold changes
  cat("Shrinking log2 fold changes...\n")
  coef_names <- grep("cluster", resultsNames(dds), value = TRUE)
  if (length(coef_names) == 0) {
    stop("No 'cluster' coefficient found for LFC shrinkage")
  }
  coef_name <- coef_names[1]
  if (length(coef_names) > 1) {
    message("Multiple 'cluster' coefficients found; using: ", coef_name)
  }
  res_shrink <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")
  res_shrink <- as.data.frame(res_shrink)

  res_shrink$sig <- rep('Not significant', nrow(res_shrink))
  res_shrink$sig[which(res_shrink$padj < alpha)] <- 'Significant'

  ## Attach pct_in using CalcPctInOutByCells
  cells_in <- colnames(seurat)[seurat@meta.data[[clus_ident]] == group_1]
  cells_out <- colnames(seurat)[seurat@meta.data[[clus_ident]] == group_2]

  pct_df <- CalcPctInOutByCells(
    seurat = seurat,
    cells_in = cells_in,
    cells_out = cells_out,
    assay = assay,
    slot = "counts"
  )

  res_shrink$feature <- rownames(res_shrink)

  merged_res <- merge(
    x = res_shrink,
    y = pct_df[, c("feature", "pct_in", "pct_out")],
    by = "feature",
    all.x = TRUE,
    sort = FALSE
  )

  res_shrink <- merged_res

  # Set up data for top gene plot
  top_gene_idx <- which.min(res$padj)
  if (length(top_gene_idx) > 0 && !is.na(top_gene_idx)) {
    d <- plotCounts(dds, gene = top_gene_idx, intgroup = "cluster", returnData = TRUE)
    exp_gene_logfc <- res[top_gene_idx, ]$log2FoldChange
    exp_gene_padj <- res[top_gene_idx, ]$padj
  } else {
    d <- NULL
    exp_gene_logfc <- NA
    exp_gene_padj <- NA
  }

  # Plots
  cat("Generating diagnostic plots...\n")
  pdf(file.path(out_dir, "diagnostic_plots.pdf"), width = 12, height = 10)

  # Pseudobulk counts distribution
  gg_counts <- cluster_counts[, sort(colnames(cluster_counts)), drop = FALSE]
  gg_counts <- tidyr::pivot_longer(
    log10(gg_counts + 1),
    cols = everything(),
    names_to = "key",
    values_to = "value"
  )
  gg_counts$condition <- cluster_metadata[gg_counts$key, "cluster"]

  print(
    ggplot(gg_counts, aes(x = key, y = value, fill = condition)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "bottom") +
      ylab("Log10(Counts + 1)") +
      ggtitle(paste(group_1, "vs", group_2, "- Pseudobulk counts"))
  )

  # PCA by condition
  print(
    DESeq2::plotPCA(vst, intgroup = "cluster") +
      theme_classic() +
      ggtitle(paste(group_1, "vs", group_2, "- PCA by condition"))
  )

  # PCA by sample
  print(
    DESeq2::plotPCA(vst, intgroup = "sample") +
      theme_classic() +
      geom_text_repel(aes(label = sample), show.legend = FALSE) +
      theme(legend.position = "bottom",
            plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      ggtitle(paste(group_1, "vs", group_2, "- PCA by sample"))
  )

  # PCA by batch if available
  if (!is.null(batch_var)) {
    print(
      DESeq2::plotPCA(vst, intgroup = batch_var) +
        theme_classic() +
        ggtitle(paste(group_1, "vs", group_2, "- PCA by batch"))
    )
  }

  # Dispersion estimates
  plotDispEsts(dds, main = paste(group_1, "vs", group_2, "- Dispersion estimates"))

  # MA plot
  plotMA(res, main = paste(group_1, "vs", group_2, "- MA plot"))

  # Volcano
  print(
    ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("grey40", "blue")) +
      theme_classic() +
      ggtitle(paste(group_1, "vs", group_2, "- Volcano plot")) +
      geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(-1, 1),      linetype = "dashed", color = "grey")
  )

  # Top gene expression
  if (!is.null(d)) {
    print(
      ggplot(d, aes(x = cluster, y = count, color = cluster)) +
        geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
        theme_classic() +
        labs(
          title    = paste("Top gene:", row.names(res)[top_gene_idx]),
          subtitle = paste0(
            "Log2FC = ", round(exp_gene_logfc, 3),
            " (", group_1, "/", group_2, ")",
            ", padj = ", formatC(exp_gene_padj, format = "e", digits = 2)
          )
        ) +
        ylab("Normalized counts") +
        theme(legend.position = "none")
    )
  }

  dev.off()

  # Save results
  write.csv(
    res_shrink,
    file = file.path(out_dir, "results.csv"),
    row.names = FALSE
  )

  # Save significant results
  res_sig <- res_shrink[which(res_shrink$padj < alpha), ]
  if (nrow(res_sig) > 0) {
    res_sig <- res_sig[order(res_sig$padj), ]
    write.csv(
      res_sig,
      file = file.path(out_dir, "significant_results.csv"),
      row.names = FALSE
    )
    cat(paste("Saved", nrow(res_sig), "significant genes\n"))
  }

  # Summary
  end <- Sys.time()
  cat("\n=== Analysis Complete ===\n")
  cat(paste("Start time:", start, "\n"))
  cat(paste("End time:", end, "\n"))
  cat(paste("Duration:", round(difftime(end, start, units = "mins"), 2), "minutes\n"))
  cat(paste("Comparison direction:", output_name, "\n"))
  cat(paste("(Positive log2FC = higher in", group_2, ")\n"))

  # Return results
  invisible(list(
    results = res_shrink,
    significant_results = if(nrow(res_sig) > 0) res_sig else NULL,
    design_formula = design_formula,
    reduced_formula = if(test_type == "LRT") reduced_formula else NULL,
    groups = c(group_1, group_2),
    dds = dds
  ))
}
