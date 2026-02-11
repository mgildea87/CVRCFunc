#' Find all markers via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param expfilt_freq genes that have greater than \code{expfilt_counts} in greater than \code{expfilt_freq} fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @param expfilt_counts genes with less than \code{expfilt_counts} in \code{expfilt_freq * sample number} will be removed from DESeq2 mode. 1 by default.
#' @param n_top_genes number of top genes per cluster to save and make a heatmap with
#' @param pct.in Filter threshold for top marker genes per cluster. For a given gene and cluster, if the fraction of cells with counts is less than pct.in it is removed from top_markers.
#' @param out_dir Name of output directory
#' @param alpha FDR adjusted p-value threshold for significance in plotting. 0.1 by default.
#' @param assay Which assay to use. RNA by default. I added this parameter to enable use of ADT data when desired.
#' @return .csv files with marker genes per \code{clus_ident}. .pdf files with plots
#' @import Seurat pheatmap DESeq2 reshape2 ggplot2 ggrepel stringr utils grDevices grr
#' @importFrom BiocGenerics t
#' @importFrom presto wilcoxauc
#' @export

FindMarkersBulk <- function(seurat,
                            clus_ident,
                            sample_ident,
                            batch_var = NULL,
                            covariates = NULL,
                            design_formula = NULL,
                            expfilt_counts = 1,
                            expfilt_freq = 0.5,
                            n_top_genes = 50,
                            pct.in = 25,
                            out_dir = "FindMarkersBulk_outs",
                            alpha = 0.1,
                            assay = 'RNA',
                            test_type = "LRT",  # "LRT" or "Wald"
                            contrast_level = NULL){  # For Wald test

  start <- Sys.time()
  coef <- variable <- value <- NULL

  # Create output directory
  dir.create(out_dir, showWarnings = FALSE)

  # Set identities
  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  # Validate inputs
  if (!is.null(batch_var) && !batch_var %in% colnames(seurat@meta.data)) {
    stop(paste("batch_var", batch_var, "not found in seurat metadata"))
  }

  if (!is.null(covariates)) {
    missing_covs <- covariates[!covariates %in% colnames(seurat@meta.data)]
    if (length(missing_covs) > 0) {
      stop(paste("Covariates not found in metadata:", paste(missing_covs, collapse = ", ")))
    }
  }

  # Build design formula if not provided
  if (is.null(design_formula)) {
    design_terms <- "iscluster"

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
      if (length(design_terms) > 1) {
        reduced_formula <- as.formula(paste("~", paste(design_terms[-1], collapse = " + ")))
      } else {
        reduced_formula <- as.formula("~ 1")
      }
    }
  } else {
    # User provided custom formula
    cat("Using custom design formula\n")
    if (test_type == "LRT") {
      # Try to create reduced formula by removing iscluster
      reduced_formula <- update(design_formula, ~ . - iscluster)
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
  pdf(paste0(out_dir, '/cells_per_clus_HM.pdf'))
  pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, sample_ident]),
           display_numbers = TRUE,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize_number = 4,
           main = "Cells per cluster and sample")
  dev.off()

  # If batch variable provided, show batch distribution
  if (!is.null(batch_var)) {
    pdf(paste0(out_dir, '/cells_per_batch_HM.pdf'))
    pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, batch_var]),
             display_numbers = TRUE,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             fontsize_number = 4,
             main = "Cells per cluster and batch")
    dev.off()
  }

  # Run wilcoxauc from presto for pct.in and pct.out
  cat("Running Wilcoxon test for pct calculations...\n")
  wilcox <- wilcoxauc(seurat, assay = 'data', seurat_assay = assay, group_by = clus_ident)

  # Initialize storage for top markers
  top_markers <- vector()
  all_results <- list()

  # Loop through each cluster
  for (cluster in clusters) {
    cat(paste("\n=== Processing cluster:", cluster, "===\n"))

    # Create cluster identifier
    iscluster <- rep(cluster, ncol(seurat))
    iscluster[which(seurat[[clus_ident]] != cluster)] <- 'other'

    # Create explicit cluster_sample identifier
    seurat$cluster_sample <- paste(iscluster,
                                   seurat@meta.data[[sample_ident]],
                                   sep = "_")

    # Pseudobulk aggregation
    cat("Aggregating counts to pseudobulk...\n")
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

    # Build metadata for pseudobulk samples
    cluster_metadata <- data.frame(
      sample_ident = sapply(stringr::str_split(colnames(cluster_counts),
                                               pattern = "_", n = 2), `[`, 2),
      iscluster = sapply(stringr::str_split(colnames(cluster_counts),
                                            pattern = "_", n = 2), `[`, 1)
    )

    # Add batch and covariate information from original metadata
    variables_to_add <- c(batch_var, covariates)
    if (!is.null(variables_to_add)) {
      # Create lookup table from original seurat object
      lookup_vars <- c(sample_ident, variables_to_add)
      sample_lookup <- unique(seurat@meta.data[, lookup_vars, drop = FALSE])

      # Merge with cluster metadata
      cluster_metadata <- merge(cluster_metadata,
                                sample_lookup,
                                by.x = "sample_ident",
                                by.y = sample_ident,
                                sort = FALSE)

      # Ensure proper ordering
      cluster_metadata <- cluster_metadata[match(colnames(cluster_counts),
                                                 paste(cluster_metadata$iscluster,
                                                       cluster_metadata$sample_ident,
                                                       sep = "_")), ]
    }

    # Set iscluster as factor with 'other' as reference
    cluster_metadata$iscluster <- factor(cluster_metadata$iscluster,
                                         levels = c('other', as.character(cluster)))

    # Convert character covariates to factors
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

    cat("Pseudobulk samples:", nrow(cluster_metadata), "\n")
    cat("  - Cluster samples:", sum(cluster_metadata$iscluster == cluster), "\n")
    cat("  - Other samples:", sum(cluster_metadata$iscluster == "other"), "\n")

    # Check for confounding
    if (!is.null(batch_var)) {
      batch_table <- table(cluster_metadata$iscluster, cluster_metadata[[batch_var]])
      cat("\nBatch distribution:\n")
      print(batch_table)

      # Warn if batches are confounded with cluster
      if (any(rowSums(batch_table > 0) == 1)) {
        warning(paste("Cluster", cluster, "may be confounded with batch variable!"))
      }
    }

    # Create DESeq2 dataset
    tryCatch({
      dds <- DESeqDataSetFromMatrix(cluster_counts,
                                    colData = cluster_metadata,
                                    design = design_formula)

      # Filter low-expressed genes
      cat(paste("Filtering genes: min", expfilt_counts, "counts in",
                expfilt_freq * 100, "% of samples\n"))
      keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq * nrow(cluster_metadata)
      cat(paste("Genes before filtering:", nrow(dds), "\n"))
      dds <- dds[keep, ]
      cat(paste("Genes after filtering:", nrow(dds), "\n"))

      if (nrow(dds) < 10) {
        warning(paste("Very few genes remaining for cluster", cluster, "- skipping"))
        next
      }

      # Variance stabilizing transformation
      cat("Performing variance stabilizing transformation...\n")
      vst <- varianceStabilizingTransformation(dds, blind = FALSE)

      # Run DESeq2
      cat(paste("Running DESeq2 with", test_type, "test...\n"))
      if (test_type == "LRT") {
        dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
        res <- results(dds, alpha = alpha)
      } else if (test_type == "Wald") {
        dds <- DESeq(dds, test = "Wald")
        # Extract results for iscluster coefficient
        contrast_name <- paste0("iscluster", cluster, "_vs_other")
        res <- results(dds,
                       name = resultsNames(dds)[grep("iscluster", resultsNames(dds))[1]],
                       alpha = alpha)
      }

      cat(paste("Significant genes (padj <", alpha, "):", sum(res$padj < alpha, na.rm = TRUE), "\n"))

      # Shrink log2 fold changes
      cat("Shrinking log2 fold changes...\n")
      coef_name <- grep("iscluster", resultsNames(dds), value = TRUE)[1]
      res_shrink <- lfcShrink(dds,
                              coef = coef_name,
                              res = res,
                              type = "apeglm")

      res_shrink <- as.data.frame(res_shrink)
      res_shrink$sig <- rep('Not significant', nrow(res_shrink))
      res_shrink$sig[which(res_shrink$padj < alpha)] <- 'Significant'

      # Set up data for top gene plot
      top_gene <- rownames(res)[which.min(res$padj)]
      if (length(top_gene) > 0 && !is.na(top_gene)) {
        d <- plotCounts(dds, gene = top_gene, intgroup = "iscluster", returnData = TRUE)
        exp_gene_logfc <- res[top_gene, ]$log2FoldChange
        exp_gene_padj <- res[top_gene, ]$padj
      } else {
        d <- NULL
        exp_gene_logfc <- NA
        exp_gene_padj <- NA
      }

      # Generate diagnostic plots
      cat("Generating diagnostic plots...\n")
      gg_counts <- cluster_counts[, sort(colnames(cluster_counts))]
      gg_counts <- melt(log10(gg_counts + 1))

      pdf(paste0(out_dir, '/cluster_', cluster, "_diagnostic_plots.pdf"))

      # Boxplot of pseudobulk counts
      print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) +
              geom_boxplot() +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1),
                    legend.position = "none") +
              ylab('Log10(Counts + 1)') +
              ggtitle(paste("Cluster", cluster, "- Pseudobulk counts")))

      # PCA by cluster
      print(DESeq2::plotPCA(vst, intgroup = "iscluster") +
              theme_classic() +
              ggtitle(paste("Cluster", cluster, "- PCA by cluster status")))

      # PCA by sample
      print(DESeq2::plotPCA(vst, intgroup = "sample_ident") +
              theme_classic() +
              geom_text_repel(aes(label = sample_ident), show.legend = FALSE) +
              ggtitle(paste("Cluster", cluster, "- PCA by sample")))

      # PCA by batch if available
      if (!is.null(batch_var)) {
        print(DESeq2::plotPCA(vst, intgroup = batch_var) +
                theme_classic() +
                ggtitle(paste("Cluster", cluster, "- PCA by batch")))
      }

      # Dispersion estimates
      plotDispEsts(dds, main = paste("Cluster", cluster, "- Dispersion estimates"))

      # MA plot
      plotMA(res, main = paste("Cluster", cluster, "- MA plot"))

      # Volcano plot
      print(ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
              geom_point(alpha = 0.7) +
              scale_color_manual(values = c('grey40', 'blue')) +
              theme_classic() +
              ggtitle(paste("Cluster", cluster, "- Volcano plot")) +
              geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
              geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey"))

      # Top gene expression
      if (!is.null(d)) {
        print(ggplot(d, aes(x = iscluster, y = count, color = iscluster)) +
                geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
                theme_classic() +
                ggtitle(paste('Top gene:', top_gene,
                              '\nLog2FC =', round(exp_gene_logfc, 3),
                              '\npadj =', formatC(exp_gene_padj, format = "e", digits = 2))) +
                ylab("Normalized counts") +
                scale_color_manual(values = c("grey60", "red")))
      }

      dev.off()

      # Add pct.in and pct.out from Wilcoxon test
      wilcox_sub <- wilcox[which(wilcox$group == cluster), ]
      res_shrink$feature <- row.names(res_shrink)
      merged_res <- merge(x = res_shrink, y = wilcox_sub, by = 'feature', sort = FALSE)
      res_shrink$pct_in <- merged_res$pct_in
      res_shrink$pct_out <- merged_res$pct_out

      # Save full results table
      write.csv(file = paste0(out_dir, '/cluster_', cluster, "_results.csv"),
                res_shrink, row.names = FALSE)

      # Store results
      all_results[[as.character(cluster)]] <- res_shrink

      # Extract top significant markers
      res_sig <- res_shrink[which(!is.na(res_shrink$padj) &
                                    res_shrink$padj < alpha &
                                    res_shrink$pct_in > pct.in), ]

      if (nrow(res_sig) > 0) {
        res_sig_sorted <- res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ]
        n_genes_to_add <- min(n_top_genes, nrow(res_sig_sorted))
        top_markers <- c(top_markers, res_sig_sorted$feature[1:n_genes_to_add])
        cat(paste("Added", n_genes_to_add, "top markers\n"))
      } else {
        cat("No significant markers found for this cluster\n")
      }

    }, error = function(e) {
      cat(paste("Error processing cluster", cluster, ":", e$message, "\n"))
      warning(paste("Skipping cluster", cluster, "due to error"))
    })
  }

  # Save all top markers
  top_markers <- unique(top_markers[!is.na(top_markers)])
  write.csv(top_markers,
            file = paste0(out_dir, '/Top_markers.csv'),
            row.names = FALSE,
            quote = FALSE)

  cat(paste("\nTotal unique top markers:", length(top_markers), "\n"))

  # Generate heatmap of top markers if any found
  if (length(top_markers) > 0) {
    cat("Generating heatmap of top markers...\n")

    # Calculate average expression
    avg_exp <- AverageExpression(seurat,
                                 assays = assay,
                                 slot = 'data',
                                 group.by = clus_ident,
                                 features = top_markers)
    avg_exp <- avg_exp[[assay]]

    # Create heatmap
    paletteLength <- 50
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    myBreaks <- c(seq(-2, 0, length.out = ceiling(paletteLength/2) + 1),
                  seq(2/paletteLength, 2, length.out = floor(paletteLength/2)))

    pdf(file = paste0(out_dir, '/Top_markers_HM.pdf'), height = 25, width = 10)
    print(pheatmap(avg_exp,
                   cluster_rows = FALSE,
                   scale = 'row',
                   color = myColor,
                   breaks = myBreaks,
                   border_color = NA,
                   main = 'Scaled average expression of top markers',
                   fontsize = 6))
    dev.off()
  } else {
    cat("No top markers to plot\n")
  }

  # Print timing
  end <- Sys.time()
  cat("\n=== Analysis Complete ===\n")
  cat(paste("Start time:", start, "\n"))
  cat(paste("End time:", end, "\n"))
  cat(paste("Duration:", round(difftime(end, start, units = "mins"), 2), "minutes\n"))

  # Return results invisibly
  invisible(list(
    all_results = all_results,
    top_markers = top_markers,
    design_formula = design_formula,
    reduced_formula = if(test_type == "LRT") reduced_formula else NULL
  ))
}
