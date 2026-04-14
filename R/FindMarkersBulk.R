#' Find all markers via pseudobulking and DESeq2
#'
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param batch_var Optional batch variable in metadata
#' @param covariates Optional vector of additional covariate names in metadata
#' @param design_formula Optional custom design formula. If NULL, will construct ~ iscluster + batch_var + covariates
#' @param expfilt_counts genes with at least \code{expfilt_counts} in
#'        \code{expfilt_freq * n_samples} will be kept. Default: 1
#' @param expfilt_freq fraction of samples that must pass \code{expfilt_counts} threshold. Default: 0.5
#' @param n_top_genes number of top genes per cluster to save and make a heatmap with
#' @param pct.in Filter threshold for top marker genes per cluster.
#'        Interpreted as percent if > 1 (e.g. 25 means 25\%). If <= 1, treated as fraction. Default: 25. Output is a fraction.
#' @param out_dir Name of output directory
#' @param alpha FDR adjusted p-value threshold for significance in plotting. Default: 0.1
#' @param assay Which assay to use. Default: 'RNA'
#' @param test_type "LRT" (likelihood ratio test) or "Wald"
#' @param contrast_level unused currently, reserved for future contrast control
#'
#' @return Invisibly returns a list with:
#'   \item{all_results}{A list of data.frames with DE results per cluster}
#'   \item{top_markers}{Vector of unique top marker genes used in heatmap}
#'   \item{design_formula}{The design formula used}
#'   \item{reduced_formula}{Reduced formula used for LRT (if applicable)}
#'   \item{params}{List of key parameters}
#'
#' outputs:
#'   - .csv files with marker genes per cluster
#'   - .pdf files with QC plots and top marker heatmap
#'
#' @import Seurat pheatmap DESeq2 ggplot2 ggrepel stringr utils grDevices grr
#' @importFrom BiocGenerics t
#' @importFrom tidyr pivot_longer
#' @export

FindMarkersBulk <- function(seurat,
                            clus_ident,
                            sample_ident,
                            batch_var = NULL,
                            covariates = NULL,
                            design_formula = NULL,
                            expfilt_counts = 1,
                            expfilt_freq = 0.5,
                            n_top_genes = 15,
                            pct.in = 25,
                            out_dir = "FindMarkersBulk_outs",
                            alpha = 0.1,
                            assay = "RNA",
                            test_type = "LRT",
                            contrast_level = NULL) {
  start <- Sys.time()
  coef <- variable <- value <- NULL

  ## ---------------------------
  ## 0. Basic input validation
  ## ---------------------------
  if (!inherits(seurat, "Seurat")) {
    stop("'seurat' must be a Seurat object")
  }
  if (!clus_ident %in% colnames(seurat@meta.data)) {
    stop("'clus_ident' not found in seurat@meta.data")
  }
  if (!sample_ident %in% colnames(seurat@meta.data)) {
    stop("'sample_ident' not found in seurat@meta.data")
  }
  if (!assay %in% names(seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object"))
  }
  if (!test_type %in% c("LRT", "Wald")) {
    stop("test_type must be one of 'LRT' or 'Wald'")
  }

  # Validate batch_var & covariates
  if (!is.null(batch_var) && !batch_var %in% colnames(seurat@meta.data)) {
    stop(paste("batch_var", batch_var, "not found in seurat metadata"))
  }
  if (!is.null(covariates)) {
    missing_covs <- covariates[!covariates %in% colnames(seurat@meta.data)]
    if (length(missing_covs) > 0) {
      stop(paste("Covariates not found in metadata:",
                 paste(missing_covs, collapse = ", ")))
    }
  }

  # set default assay
  DefaultAssay(seurat) <- assay

  #Normalize
  seurat <- NormalizeData(seurat, assay = assay)

  # Interpret pct.in as percent or fraction
  pct_in_threshold <- if (pct.in > 1) pct.in / 100 else pct.in

  # Create output directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Set identities
  Idents(seurat) <- clus_ident
  clusters <- sort(as.vector(unique(Idents(seurat))))

  ## ---------------------------
  ## 1. Design formula helpers
  ## ---------------------------
  .build_design <- function(batch_var, covariates) {
    terms <- "iscluster"
    if (!is.null(batch_var)) {
      terms <- c(terms, batch_var)
      message("Adding batch variable to design: ", batch_var)
    }
    if (!is.null(covariates)) {
      terms <- c(terms, covariates)
      message("Adding covariates to design: ", paste(covariates, collapse = ", "))
    }
    as.formula(paste("~", paste(terms, collapse = " + ")))
  }

  .build_reduced <- function(design_formula) {
    # Remove iscluster, if nothing left then ~1
    reduced <- update(design_formula, ~ . - iscluster)
    if (length(attr(terms(reduced), "term.labels")) == 0) {
      reduced <- ~ 1
    }
    reduced
  }

  # Build or validate design formula
  if (is.null(design_formula)) {
    design_formula <- .build_design(batch_var, covariates)
    if (test_type == "LRT") {
      reduced_formula <- .build_reduced(design_formula)
    } else {
      reduced_formula <- NULL
    }
  } else {
    message("Using custom design formula")
    design_terms <- attr(terms(design_formula), "term.labels")
    if (!any(grepl("^iscluster", design_terms))) {
      stop("Design formula must include 'iscluster' (e.g. ~ iscluster + batch)")
    }
    if (test_type == "LRT") {
      reduced_formula <- .build_reduced(design_formula)
    } else {
      reduced_formula <- NULL
    }
  }

  message("\nFull design formula: ", deparse(design_formula))
  if (test_type == "LRT") {
    message("Reduced design formula: ", deparse(reduced_formula))
  }
  message("Test type: ", test_type, "\n")

  ## ---------------------------
  ## 2. QC heatmaps of cell numbers
  ## ---------------------------
  pdf(file.path(out_dir, "cells_per_clus_HM.pdf"))
  pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, sample_ident]),
           display_numbers = TRUE,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           fontsize_number = 4,
           main = "Cells per cluster and sample")
  dev.off()

  if (!is.null(batch_var)) {
    pdf(file.path(out_dir, "cells_per_batch_HM.pdf"))
    pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, batch_var]),
             display_numbers = TRUE,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             fontsize_number = 4,
             main = "Cells per cluster and batch")
    dev.off()
  }

  ## ---------------------------
  ## 3. Initialize outputs
  ## ---------------------------
  top_markers <- character(0)
  all_results <- list()

  ## ---------------------------
  ## 4. Loop over clusters
  ## ---------------------------
  for (cluster in clusters) {
    message("\n=== Processing cluster: ", cluster, " ===")

    cluster_success <- tryCatch({
      ## 5.1 Construct iscluster and cluster_sample
      iscluster <- rep(cluster, ncol(seurat))
      iscluster[seurat[[clus_ident]][, 1] != cluster] <- "other"

      seurat$cluster_sample <- paste(iscluster, seurat[[sample_ident]][, 1], sep = "_")

      ## 5.2 Pseudobulk aggregation (version-agnostic)
      message("Aggregating counts to pseudobulk...")
      counts_mat <- GetAssayData(seurat, assay = assay, slot = "counts")
      #if (!is.matrix(counts_mat)) {
      #  counts_mat <- as.matrix(counts_mat)
      #}
      groupings <- data.frame(cluster_sample = seurat$cluster_sample)
      pb <- aggregate.Matrix(
        t(counts_mat),
        groupings = groupings,
        fun = "sum"
      )

      cluster_sample_names <- row.names(pb)

      ## 5.3 Build cluster_metadata
      cluster_metadata <- data.frame(
        iscluster = sub("_.*", "", cluster_sample_names),
        sample_ident = sub("^[^_]*_", "", cluster_sample_names),
        stringsAsFactors = FALSE
      )

      vars_to_add <- unique(c(sample_ident, batch_var, covariates))
      sample_meta <- unique(seurat@meta.data[, vars_to_add, drop = FALSE])
      colnames(sample_meta)[colnames(sample_meta) == sample_ident] <- "sample_ident"

      cluster_metadata <- merge(
        cluster_metadata,
        sample_meta,
        by = "sample_ident",
        all.x = TRUE,
        sort = FALSE
      )

      # Ensure row order matches columns of pb
      cluster_metadata <- cluster_metadata[match(
        cluster_sample_names,
        paste(cluster_metadata$iscluster, cluster_metadata$sample_ident, sep = "_")
      ), ]
      rownames(cluster_metadata) <- cluster_sample_names

      ## 5.4 Factor conversions
      cluster_metadata$iscluster <- factor(
        cluster_metadata$iscluster,
        levels = c("other", as.character(cluster))
      )

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

      message("Pseudobulk samples: ", nrow(cluster_metadata))
      message("  - Cluster samples: ",
              sum(cluster_metadata$iscluster == cluster, na.rm = TRUE))
      message("  - Other samples: ",
              sum(cluster_metadata$iscluster == "other", na.rm = TRUE))

      # Check minimal samples for each group
      if (sum(cluster_metadata$iscluster == cluster) < 2 ||
          sum(cluster_metadata$iscluster == "other") < 2) {
        warning("Too few samples in cluster or other for cluster ", cluster, " - skipping.")
        return(FALSE)
      }

      ## 5.5 Check for confounding batch
      if (!is.null(batch_var)) {
        batch_table <- table(cluster_metadata$iscluster, cluster_metadata[[batch_var]])
        message("\nBatch distribution:")
        print(batch_table)
        if (any(rowSums(batch_table > 0) == 1)) {
          warning(paste("Cluster", cluster, "may be confounded with batch variable!"))
        }
      }

      ## 5.6 Create DESeq2 object
      cluster_counts <- as.data.frame(t(as.matrix(pb)))

      dds <- DESeqDataSetFromMatrix(
        countData = cluster_counts,
        colData   = cluster_metadata,
        design    = design_formula
      )

      ## 5.7 Gene filtering
      message("Filtering genes: min ", expfilt_counts, " counts in ",
              expfilt_freq * 100, "% of samples")
      keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq * nrow(cluster_metadata)
      message("Genes before filtering: ", nrow(dds))
      dds <- dds[keep, ]
      message("Genes after filtering: ", nrow(dds))

      if (nrow(dds) < 10) {
        warning(paste("Very few genes remaining for cluster", cluster, "- skipping"))
        return(FALSE)
      }

      ## 5.8 VST and DESeq2
      message("Performing variance stabilizing transformation...")
      vst <- varianceStabilizingTransformation(dds, blind = FALSE)

      message("Running DESeq2 with ", test_type, " test...")
      if (test_type == "LRT") {
        dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
        res <- results(dds, alpha = alpha)
      } else { # Wald
        dds <- DESeq(dds, test = "Wald")
        iscluster_coef <- grep("iscluster", resultsNames(dds))[1]
        if (is.na(iscluster_coef)) {
          stop("No 'iscluster' coefficient found in resultsNames(dds)")
        }
        res <- results(
          dds,
          name  = resultsNames(dds)[iscluster_coef],
          alpha = alpha
        )
      }

      message("Significant genes (padj < ", alpha, "): ",
              sum(res$padj < alpha, na.rm = TRUE))

      ## 5.9 LFC shrinkage with fallback
      message("Shrinking log2 fold changes...")
      coef_names <- grep("iscluster", resultsNames(dds), value = TRUE)
      if (length(coef_names) == 0) {
        stop("No 'iscluster' coefficient found for LFC shrinkage in cluster ", cluster)
      }
      coef_name <- coef_names[1]
      if (length(coef_names) > 1) {
        message("Multiple 'iscluster' coefficients found; using: ", coef_name)
      }

      res_shrink <- tryCatch(
        {
          lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")
        },
        error = function(e) {
          message("apeglm shrinkage failed for cluster ", cluster, ": ", e$message,
                  ". Falling back to 'normal' shrinkage.")
          lfcShrink(dds, coef = coef_name, res = res, type = "normal")
        }
      )
      res_shrink <- as.data.frame(res_shrink)
      res_shrink$sig <- "Not significant"
      res_shrink$sig[which(res_shrink$padj < alpha)] <- "Significant"

      ## 5.10 Top gene
      res <- res[order(res$log2FoldChange, decreasing = T),]
      res_sig <- res[which(res$padj < 0.05),]
      top_gene <- rownames(res_sig)[1]
      if (!is.null(top_gene) && length(top_gene) > 0 && !is.na(top_gene)) {
        d <- plotCounts(dds, gene = top_gene, intgroup = "iscluster", returnData = TRUE)
        exp_gene_logfc <- res[top_gene, ]$log2FoldChange
        exp_gene_padj  <- res[top_gene, ]$padj
      } else {
        d <- NULL
        exp_gene_logfc <- NA
        exp_gene_padj  <- NA
      }

      ## 5.11 Diagnostic plots
      message("Generating diagnostic plots...")
      pdf_file <- file.path(out_dir, paste0("cluster_", cluster, "_diagnostic_plots.pdf"))
      pdf(pdf_file)

      # pseudobulk counts distribution
      gg_counts <- cluster_counts[, sort(colnames(cluster_counts)), drop = FALSE]
      gg_counts <- tidyr::pivot_longer(
        log10(gg_counts + 1),
        cols = everything(),
        names_to = "key",
        values_to = "value"
      )
      gg_counts$type <- sub("_.*", "", gg_counts$key)
      gg_counts$sample <- sub("^[^_]*_", "", gg_counts$key)

      print(
        ggplot(gg_counts, aes(x = sample, y = value, fill = type)) +
          geom_boxplot() +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = "none") +
          ylab("Log10(Counts + 1)") +
          ggtitle(paste("Cluster", cluster, "- Pseudobulk counts"))
      )

      # PCA by cluster
      print(
        DESeq2::plotPCA(vst, intgroup = "iscluster") +
          theme_classic() +
          ggtitle(paste("Cluster", cluster, "- PCA by cluster status"))
      )

      # PCA by sample
      print(
        DESeq2::plotPCA(vst, intgroup = "sample_ident") +
          theme_classic() +
          theme(legend.position = "bottom",
                plot.margin = unit(c(0, 0, 0, 0), "cm")) +
          ggtitle(paste("Cluster", cluster, "- PCA by sample"))
      )

      # PCA by batch if available
      if (!is.null(batch_var)) {
        print(
          DESeq2::plotPCA(vst, intgroup = batch_var) +
            theme_classic() +
            ggtitle(paste("Cluster", cluster, "- PCA by batch"))
        )
      }

      # Dispersion estimates
      plotDispEsts(dds, main = paste("Cluster", cluster, "- Dispersion estimates"))

      # MA plot
      plotMA(res, main = paste("Cluster", cluster, "- MA plot"))

      # Volcano
      print(
        ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
          geom_point(alpha = 0.7) +
          scale_color_manual(values = c("grey40", "blue")) +
          theme_classic() +
          ggtitle(paste("Cluster", cluster, "- Volcano plot")) +
          geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey")
      )

      # Top gene expression
      if (!is.null(d)) {
        print(
          ggplot(d, aes(x = iscluster, y = count, color = iscluster)) +
            geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
            theme_classic() +
            labs(
              title = paste("Top gene:", top_gene),
              subtitle = paste0(
                "Log2FC = ", round(exp_gene_logfc, 3),
                ", padj = ", formatC(exp_gene_padj, format = "e", digits = 2)
              )
            ) +
            ylab("Normalized counts") +
            scale_color_manual(values = c("grey60", "red"))
        )
      }

      # Close device for this cluster
      dev.off()

      ## 5.12 Attach pct_in / pct_out
      ## Compute pct_in / pct_out for this cluster using cell barcodes
      cells_in  <- colnames(seurat)[seurat[[clus_ident]][, 1] == cluster]
      cells_out <- colnames(seurat)[seurat[[clus_ident]][, 1] != cluster]

      pct_df_cluster <- CalcPctInOutByCells(
        seurat   = seurat,
        cells_in = cells_in,
        cells_out = cells_out,
        assay    = assay,
        slot     = "counts"   # or "counts" if you prefer
      )

      # Merge into DESeq2 results
      res_shrink$feature <- rownames(res_shrink)

      merged_res <- merge(
        x   = res_shrink,
        y   = pct_df_cluster[, c("feature", "pct_in", "pct_out")],
        by  = "feature",
        all.x = TRUE,
        sort  = FALSE
      )

      res_shrink <- merged_res

      ## 5.13 Save full results
      write.csv(
        res_shrink,
        file = file.path(out_dir, paste0("cluster_", cluster, "_results.csv")),
        row.names = FALSE
      )

      all_results[[as.character(cluster)]] <- res_shrink

      ## 5.14 Select top markers for cluster
      res_sig <- subset(
        res_shrink,
        !is.na(padj) &
          padj < alpha &
          !is.na(pct_in) &
          pct_in > pct_in_threshold &
          log2FoldChange > 0
      )

      if (nrow(res_sig) > 0) {
        res_sig_sorted <- res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ]
        n_genes_to_add <- min(n_top_genes, nrow(res_sig_sorted))
        top_markers <- c(top_markers, res_sig_sorted$feature[1:n_genes_to_add])
        message("Added ", n_genes_to_add, " top markers")
      } else {
        message("No significant markers found for this cluster")
      }

      TRUE
    }, error = function(e) {
      message("\n!!! ERROR processing cluster ", cluster, " !!!")
      message("Error message: ", e$message)
      # Close any open graphics devices
      while (dev.cur() > 1) dev.off()
      FALSE
    })
    if (!cluster_success) {
      message("*** SKIPPING cluster ", cluster, " - continuing to next cluster ***")
      next
    }
    message("Successfully completed cluster ", cluster)
  }

  ## ---------------------------
  ## 5. Top markers heatmap
  ## ---------------------------
  top_markers <- unique(top_markers[!is.na(top_markers)])
  write.csv(
    top_markers,
    file = file.path(out_dir, "Top_markers.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  message("\nTotal unique top markers: ", length(top_markers))

  if (length(top_markers) > 0) {
    message("Generating heatmap of top markers...")

    # Ensure markers exist in assay
    top_markers_valid <- intersect(top_markers, rownames(seurat[[assay]]))
    if (length(top_markers_valid) == 0) {
      message("No valid top markers present in assay '", assay, "' - skipping heatmap.")
    } else {
      # Handle Seurat v4/v5 AverageExpression API differences
      agg_formals <- names(formals(AverageExpression))

      if ("layer" %in% agg_formals) {
        avg_exp_list <- AverageExpression(
          seurat,
          assays   = assay,
          group.by = clus_ident,
          features = top_markers_valid,
          layer    = "data"
        )
      } else {
        avg_exp_list <- AverageExpression(
          seurat,
          assays   = assay,
          group.by = clus_ident,
          features = top_markers_valid,
          slot     = "data"
        )
      }

      # Seurat typically returns a list indexed by assay
      avg_exp <- avg_exp_list[[assay]]

      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(
        seq(-2, 0, length.out = ceiling(paletteLength / 2) + 1),
        seq(2 / paletteLength, 2, length.out = floor(paletteLength / 2))
      )

      pdf(file = file.path(out_dir, "Top_markers_HM.pdf"), height = 25, width = 10)
      print(
        pheatmap(avg_exp,
                 cluster_rows = TRUE,
                 scale = "row",
                 color = myColor,
                 breaks = myBreaks,
                 border_color = NA,
                 main = "Scaled average expression of top markers",
                 fontsize = 6)
      )
      dev.off()
    }
  } else {
    message("No top markers to plot")
  }

  ## ---------------------------
  ## 6. Timing and return
  ## ---------------------------
  end <- Sys.time()
  message("\n=== Analysis Complete ===")
  message("Start time: ", start)
  message("End time:   ", end)
  message("Duration:   ", round(difftime(end, start, units = "mins"), 2), " minutes")

  invisible(list(
    all_results     = all_results,
    top_markers     = top_markers,
    design_formula  = design_formula,
    reduced_formula = if (test_type == "LRT") reduced_formula else NULL,
    params          = list(
      alpha          = alpha,
      expfilt_counts = expfilt_counts,
      expfilt_freq   = expfilt_freq,
      n_top_genes    = n_top_genes,
      pct.in         = pct.in,
      assay          = assay,
      test_type      = test_type,
      batch_var      = batch_var,
      covariates     = covariates
    )
  ))
}
