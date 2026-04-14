#' Differential expression testing between conditions via pseudobulking and DESeq2
#'
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param condition_ident Identity class for conditions to be tested
#' @param conditions A vector of exactly 2 conditions within \code{condition_ident} to compare.
#'        The first element is the numerator (positive log2FC = higher in conditions[1]).
#'        If \code{NULL} and \code{test_type = "LRT"}, all levels of \code{condition_ident} are used.
#' @param batch_var Optional batch variable in metadata
#' @param covariates Optional vector of additional covariate names in metadata
#' @param design_formula Optional custom design formula. If NULL, will construct ~ condition_ident + batch_var + covariates. If providing a custom formula, all variables (except condition_ident) must be specified as batch_var or covariates to ensure they are included in the pseudobulk metadata.
#' @param reduced_formula Optional custom reduced formula for LRT. If NULL, will construct by removing the condition_ident term from design_formula
#' @param expfilt_counts Genes with at least \code{expfilt_counts} counts in
#'        \code{expfilt_freq * n_samples} samples will be kept. Default: 1
#' @param expfilt_freq Fraction of samples that must pass \code{expfilt_counts} threshold. Default: 0.5
#' @param n_top_genes Number of top genes per cluster to save and make a heatmap with
#' @param pct.in Filter threshold for top marker genes. Interpreted as percent if > 1. Default: 25
#' @param out_dir Name of output directory
#' @param alpha FDR adjusted p-value threshold for significance in plotting. Default: 0.1
#' @param assay Which assay to use. Default: 'RNA'
#' @param test_type "LRT" (likelihood ratio test) or "Wald"
#' @param contrast_level Unused currently, reserved for future contrast control
#'
#' @return Invisibly returns a list with:
#'   \item{all_results}{A list of data.frames with DE results per cluster}
#'   \item{summary}{Summary data.frame of significant genes per cluster}
#'   \item{top_markers}{Vector of unique top marker genes used in heatmap}
#'   \item{design_formula}{The design formula used}
#'   \item{reduced_formula}{Reduced formula used for LRT (if applicable)}
#'   \item{conditions}{The conditions compared}
#'   \item{params}{List of key parameters}
#'
#' outputs:
#'   - .csv files with DE results per cluster
#'   - .pdf files with QC plots and top marker heatmap
#'
#' @import Seurat pheatmap DESeq2 ggplot2 ggrepel stringr utils grDevices grr
#' @importFrom BiocGenerics t
#' @importFrom tidyr pivot_longer
#' @export

FindMarkersCondition <- function(seurat,
                                 clus_ident,
                                 sample_ident,
                                 condition_ident,
                                 conditions = NULL,
                                 batch_var = NULL,
                                 covariates = NULL,
                                 design_formula = NULL,
                                 reduced_formula = NULL,
                                 expfilt_counts = 1,
                                 expfilt_freq = 0.5,
                                 n_top_genes = 15,
                                 pct.in = 25,
                                 out_dir = "FindMarkersCondition_outs",
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
  if (!condition_ident %in% colnames(seurat@meta.data)) {
    stop(paste("'condition_ident'", condition_ident, "not found in metadata"))
  }
  if (!assay %in% names(seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object"))
  }
  if (!test_type %in% c("LRT", "Wald")) {
    stop("test_type must be one of 'LRT' or 'Wald'")
  }
  # Handle conditions
  if (test_type == "LRT" && is.null(conditions)) {
    conditions <- sort(unique(seurat@meta.data[[condition_ident]]))
    message("Using all levels of condition_ident for LRT: ", paste(conditions, collapse = ", "))
  } else if (!is.null(conditions)) {
    if (length(conditions) != 2) {
      stop("Exactly 2 conditions must be provided when test_type = 'Wald' or when specifying conditions")
    }
    if (!all(conditions %in% unique(seurat@meta.data[[condition_ident]]))) {
      stop("One or both conditions not found in the data")
    }
  } else {
    stop("conditions must be specified when test_type = 'Wald'")
  }
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

  # Set default assay and normalize
  DefaultAssay(seurat) <- assay
  seurat <- NormalizeData(seurat, assay = assay)

  # Interpret pct.in as percent or fraction
  pct_in_threshold <- if (pct.in > 1) pct.in / 100 else pct.in

  # Create output directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Set identities
  Idents(seurat) <- clus_ident
  clusters <- sort(as.vector(unique(Idents(seurat))))

  message("Analyzing ", length(clusters), " clusters: ",
          paste(clusters, collapse = ", "))
  if (test_type == "LRT"){
    message("Running LRT across all levels of condition_ident")
  }
  if (test_type == "Wald") {
    message("Comparing conditions: ", conditions[1], " vs ", conditions[2],
          " (positive log2FC = higher in ", conditions[1], ")\n")
  }

  ## ---------------------------
  ## 1. Design formula helpers
  ## ---------------------------
  .build_design <- function(batch_var, covariates) {
    terms <- condition_ident
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
    reduced <- update(design_formula, as.formula(paste("~ . -", condition_ident)))
    if (length(attr(terms(reduced), "term.labels")) == 0) {
      reduced <- ~ 1
    }
    reduced
  }

  # Build or validate design formula
  if (is.null(design_formula)) {
    design_formula <- .build_design(batch_var, covariates)
    if (test_type == "LRT") {
      if (is.null(reduced_formula)) {
        reduced_formula <- .build_reduced(design_formula)
      }
    } else {
      reduced_formula <- NULL
    }
  } else {
    message("Using custom design formula")
    design_terms <- attr(terms(design_formula), "term.labels")
    if (!any(grepl(paste0("^", condition_ident), design_terms))) {
      stop("Design formula must include '", condition_ident, "' (e.g. ~ ", condition_ident, " + batch)")
    }
    if (test_type == "LRT") {
      if (is.null(reduced_formula)) {
        reduced_formula <- .build_reduced(design_formula)
      }
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
  ## 2. Sample-condition metadata
  ## ---------------------------
  vars_to_keep <- unique(c(sample_ident, condition_ident, batch_var, covariates))
  samplecon_table <- unique(seurat@meta.data[, vars_to_keep, drop = FALSE])
  # Rename for internal consistency
  colnames(samplecon_table)[colnames(samplecon_table) == sample_ident]    <- "sample"
  # Keep condition_ident as is for design formula flexibility
  condition_col <- condition_ident

  message("Sample-Condition Distribution:")
  print(table(samplecon_table[[condition_col]]))
  cat("\n")

  ## ---------------------------
  ## 3. QC heatmaps of cell numbers
  ## ---------------------------
  pdf(file.path(out_dir, "cells_per_clus_HM.pdf"))
  pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, sample_ident]),
           display_numbers = TRUE,
           cluster_rows    = FALSE,
           cluster_cols    = FALSE,
           fontsize_number = 4,
           main            = "Cells per cluster and sample")
  dev.off()

  if (!is.null(batch_var)) {
    pdf(file.path(out_dir, "cells_per_batch_HM.pdf"))
    pheatmap(table(seurat@meta.data[, clus_ident], seurat@meta.data[, batch_var]),
             display_numbers = TRUE,
             cluster_rows    = FALSE,
             cluster_cols    = FALSE,
             fontsize_number = 4,
             main            = "Cells per cluster and batch")
    dev.off()
  }

  ## ---------------------------
  ## 4. Initialize outputs
  ## ---------------------------
  top_markers      <- character(0)
  all_results      <- list()
  clusters_summary <- character(0)
  sig_up           <- integer(0)
  sig_down         <- integer(0)
  sig_LRT          <- integer(0)

  ## ---------------------------
  ## 5. Loop over clusters
  ## ---------------------------
  for (cluster in clusters) {
    message("\n=== Processing cluster: ", cluster, " ===")

    cluster_success <- tryCatch({

      ## 5.1 Create cluster_sample identifier
      seurat$cluster_sample <- paste(seurat@meta.data[[clus_ident]],
                                     seurat@meta.data[[sample_ident]],
                                     sep = "_")

      ## 5.2 Pseudobulk aggregation (version-agnostic)
      message("Aggregating counts to pseudobulk...")
      counts_mat <- GetAssayData(seurat, assay = assay, slot = "counts")
      groupings  <- data.frame(cluster_sample = seurat$cluster_sample)

      pb <- aggregate.Matrix(
        t(counts_mat),
        groupings = groupings,
        fun       = "sum"
      )

      # Keep only rows belonging to this cluster
      cluster_rows <- startsWith(rownames(pb), prefix = paste0(cluster, "_"))
      pb_cluster   <- pb[cluster_rows, , drop = FALSE]

      # Extract sample names from rownames (everything after first "_")
      sample_names <- sub("^[^_]+_", "", rownames(pb_cluster))

      # Filter for samples present in samplecon_table
      keep_samples   <- sample_names %in% samplecon_table$sample
      pb_cluster     <- pb_cluster[keep_samples, , drop = FALSE]
      sample_names   <- sample_names[keep_samples]

      # Transpose: genes x samples
      cluster_counts <- as.data.frame(t(as.matrix(pb_cluster)))
      colnames(cluster_counts) <- sample_names

      ## 5.3 Build cluster_metadata aligned to cluster_counts columns
      cluster_metadata <- samplecon_table[match(sample_names, samplecon_table$sample), ]
      rownames(cluster_metadata) <- sample_names

      if(test_type == "Wald"){
        # Set condition factor levels: conditions[2] as reference
        all_condition_levels <- unique(samplecon_table$condition)
        other_levels <- setdiff(all_condition_levels, conditions)

        cluster_metadata$condition <- factor(
          cluster_metadata$condition,
          levels = c(conditions[2], conditions[1], other_levels)  # reference = conditions[2]
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

        n_cond1 <- sum(cluster_metadata$condition == conditions[1], na.rm = TRUE)
        n_cond2 <- sum(cluster_metadata$condition == conditions[2], na.rm = TRUE)

        message("Pseudobulk samples - ", conditions[1], ": ", n_cond1,
                " | ", conditions[2], ": ", n_cond2)

        if (n_cond1 < 2 || n_cond2 < 2) {
          warning("Too few samples in one or both conditions for cluster ",
                  cluster, " - skipping.")
          return(FALSE)
        }
      }
      ## 5.4 Check for batch confounding
      if (!is.null(batch_var)) {
        batch_table <- table(cluster_metadata$condition, cluster_metadata[[batch_var]])
        message("\nBatch distribution:")
        print(batch_table)
        if (any(rowSums(batch_table > 0) == 1)) {
          warning(paste("Cluster", cluster,
                        "may be confounded with batch variable!"))
        }
      }

      ## 5.5 Create DESeq2 object
      dds <- DESeqDataSetFromMatrix(
        countData = cluster_counts,
        colData   = cluster_metadata,
        design    = design_formula
      )

      ## 5.6 Gene filtering
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

      ## 5.7 VST and DESeq2
      message("Performing variance stabilizing transformation...")
      vst <- varianceStabilizingTransformation(dds, blind = FALSE)

      message("Running DESeq2 with ", test_type, " test...")
      if (test_type == "LRT") {
        dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
        res <- results(dds, alpha = alpha)
      } else { # Wald
        dds <- DESeq(dds, test = "Wald")
        res <- results(dds, contrast = c(condition_ident, conditions[1], conditions[2]), alpha = alpha)
      }

      message("Saving DESeq2 object for cluster ", cluster, "...")
      saveRDS(
        dds,
        file = file.path(out_dir, paste0("cluster_", cluster, "_dds.rds"))
      )

      message("Significant genes (padj < ", alpha, "): ",
              sum(res$padj < alpha, na.rm = TRUE))

      ## 5.8 LFC shrinkage with fallback
      if (test_type == "Wald") {
        message("Shrinking log2 fold changes...")
        res_raw <- as.data.frame(res)
        coef_names <- grep(condition_ident, resultsNames(dds), value = TRUE)
        if (length(coef_names) == 0) {
          stop("No '", condition_ident, "' coefficient found for LFC shrinkage in cluster ", cluster)
        }
        coef_name <- coef_names[1]
        if (length(coef_names) > 1) {
          message("Multiple '", condition_ident, "' coefficients found; using: ", coef_name)
        }

        res_shrink <- tryCatch(
          {
            lfcShrink(
              dds,
              contrast = c(condition_ident, conditions[1], conditions[2]),
              res      = res,
              type     = "ashr"   # ashr supports contrast=; apeglm does not
            )
          },
          error = function(e) {
            message("ashr shrinkage failed for cluster ", cluster, ": ", e$message,
                    ". Falling back to 'normal' shrinkage.")
            lfcShrink(
              dds,
              contrast = c(condition_ident, conditions[1], conditions[2]),
              res      = res,
              type     = "normal"
            )
          }
        )
        res_shrink <- as.data.frame(res_shrink)
        if ("log2FoldChange" %in% colnames(res_raw)) {
          res_shrink$log2FoldChange_raw <- res_raw$log2FoldChange
        }
      } else {
        message("Skipping LFC shrinkage for LRT; using raw DESeq2 results.")
        res_shrink <- as.data.frame(res)
        if ("log2FoldChange" %in% colnames(res_shrink)) {
          res_shrink$log2FoldChange_raw <- res_shrink$log2FoldChange
        }
      }
      res_shrink$sig  <- "Not significant"
      res_shrink$sig[which(res_shrink$padj < alpha)] <- "Significant"

      ## 5.9 Top gene plot data
      res          <- res[order(res$padj, na.last = TRUE), ]
      top_gene     <- rownames(res)[1]
      if (!is.null(top_gene) && length(top_gene) > 0 && !is.na(top_gene)) {
        d               <- plotCounts(dds, gene = top_gene,
                                      intgroup = condition_ident, returnData = TRUE)
        exp_gene_logfc  <- res[top_gene, ]$log2FoldChange
        exp_gene_padj   <- res[top_gene, ]$padj
      } else {
        d               <- NULL
        exp_gene_logfc  <- NA
        exp_gene_padj   <- NA
      }

      ## 5.10 Diagnostic plots
      message("Generating diagnostic plots...")
      pdf_file <- file.path(out_dir,
                            paste0("cluster_", cluster, "_diagnostic_plots.pdf"))
      pdf(pdf_file, width = 12, height = 10)

      # Pseudobulk counts distribution
      gg_counts <- cluster_counts[, sort(colnames(cluster_counts)), drop = FALSE]
      gg_counts <- tidyr::pivot_longer(
        log10(gg_counts + 1),
        cols = everything(),
        names_to = "key",
        values_to = "value"
      )
      gg_counts$condition <- cluster_metadata[gg_counts$key, condition_ident]

      print(
        ggplot(gg_counts, aes(x = key, y = value, fill = condition)) +
          geom_boxplot() +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = "bottom") +
          ylab("Log10(Counts + 1)") +
          ggtitle(paste("Cluster", cluster, "- Pseudobulk counts"))
      )

      # PCA by condition
      print(
        DESeq2::plotPCA(vst, intgroup = condition_ident) +
          theme_classic() +
          ggtitle(paste("Cluster", cluster, "- PCA by", condition_ident))
      )

      # PCA by sample
      print(
        DESeq2::plotPCA(vst, intgroup = "sample") +
          theme_classic() +
          geom_text_repel(aes(label = sample), show.legend = FALSE) +
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
          geom_vline(xintercept = c(-1, 1),      linetype = "dashed", color = "grey")
      )

      # Top gene expression
      if (!is.null(d)) {
        print(
          ggplot(d, aes(x = !!sym(condition_ident), y = count, color = !!sym(condition_ident))) +
            geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +
            theme_classic() +
            labs(
              title    = paste("Top gene:", top_gene),
              subtitle = paste0(
                "Log2FC = ", round(exp_gene_logfc, 3),
                " (", conditions[1], "/", conditions[2], ")",
                ", padj = ", formatC(exp_gene_padj, format = "e", digits = 2)
              )
            ) +
            ylab("Normalized counts") +
            theme(legend.position = "none")
        )
      }

      dev.off()

      ## 5.11 Attach pct_in per condition using CalcPctInOutByCells
      cells_cond1 <- colnames(seurat)[
        seurat@meta.data[[clus_ident]] == cluster &
          seurat@meta.data[[condition_ident]] == conditions[1]
      ]
      cells_cond2 <- colnames(seurat)[
        seurat@meta.data[[clus_ident]] == cluster &
          seurat@meta.data[[condition_ident]] == conditions[2]
      ]

      pct_df_cond1 <- CalcPctInOutByCells(
        seurat    = seurat,
        cells_in  = cells_cond1,
        cells_out = cells_cond2,
        assay     = assay,
        slot      = "counts"
      )
      pct_df_cond2 <- CalcPctInOutByCells(
        seurat    = seurat,
        cells_in  = cells_cond2,
        cells_out = cells_cond1,
        assay     = assay,
        slot      = "counts"
      )

      res_shrink$feature <- rownames(res_shrink)

      merged_res <- merge(
        x     = res_shrink,
        y     = pct_df_cond1[, c("feature", "pct_in")],
        by    = "feature",
        all.x = TRUE,
        sort  = FALSE
      )
      colnames(merged_res)[colnames(merged_res) == "pct_in"] <-
        paste0("pct_in_", conditions[1])

      merged_res <- merge(
        x     = merged_res,
        y     = pct_df_cond2[, c("feature", "pct_in")],
        by    = "feature",
        all.x = TRUE,
        sort  = FALSE
      )
      colnames(merged_res)[colnames(merged_res) == "pct_in"] <-
        paste0("pct_in_", conditions[2])

      res_shrink <- merged_res

      ## 5.12 Save full results
      write.csv(
        res_shrink,
        file      = file.path(out_dir, paste0("cluster_", cluster, "_results.csv")),
        row.names = FALSE
      )

      all_results[[as.character(cluster)]] <- res_shrink

      # Update summary tallies
      clusters_summary <- c(clusters_summary, as.character(cluster))
      if (test_type == "Wald") {
      sig_up   <- c(sig_up,   sum(res_shrink$padj < alpha &
                                    res_shrink$log2FoldChange > 0, na.rm = TRUE))
      sig_down <- c(sig_down, sum(res_shrink$padj < alpha &
                                    res_shrink$log2FoldChange < 0, na.rm = TRUE))
      }

      if (test_type == "LRT") {
        sig_LRT   <- c(sig_LRT,   sum(res_shrink$padj < alpha, na.rm = TRUE))
      }

      ## 5.13 Select top markers for heatmap
      pct_col <- paste0("pct_in_", conditions[1])
      res_sig <- res_shrink[
        !is.na(res_shrink$padj) &
        res_shrink$padj < alpha &
        !is.na(res_shrink[[pct_col]]) &
        res_shrink[[pct_col]] > pct_in_threshold &
        res_shrink$log2FoldChange > 0,
      ]
      if (nrow(res_sig) > 0) {
        res_sig_sorted <- res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ]
        n_genes_to_add <- min(n_top_genes, nrow(res_sig_sorted))
        top_markers    <- c(top_markers, res_sig_sorted$feature[1:n_genes_to_add])
        message("Added ", n_genes_to_add, " top markers")
      } else {
        message("No significant markers found for this cluster")
      }

      TRUE
    }, error = function(e) {
      message("\n!!! ERROR processing cluster ", cluster, " !!!")
      message("Error message: ", e$message)
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
  ## 6. Summary results
  ## ---------------------------
  summary_results <- NULL
  if (test_type == "Wald") {
    if (length(clusters_summary) > 0) {
      summary_results <- data.frame(
        cluster  = clusters_summary,
        sig_up   = sig_up,
        sig_down = sig_down,
        stringsAsFactors = FALSE
      )
      colnames(summary_results)[2:3] <- c(
        paste0("sig_up_in_",   conditions[1]),
        paste0("sig_down_in_", conditions[1])
      )
      write.csv(
        summary_results,
        file      = file.path(out_dir, "Summary_results.csv"),
        quote     = FALSE,
        row.names = FALSE
      )
      message("\n=== Summary Results ===")
      print(summary_results)
    }
  }
  if (test_type == "LRT") {
    if (length(clusters_summary) > 0) {
      summary_results <- data.frame(
        cluster  = clusters_summary,
        sig_LRT   = sig_LRT,
        stringsAsFactors = FALSE
      )
      colnames(summary_results)[2] <- 'sig_LRT'
      write.csv(
        summary_results,
        file      = file.path(out_dir, "Summary_results.csv"),
        quote     = FALSE,
        row.names = FALSE
      )
      message("\n=== Summary Results ===")
      print(summary_results)
    }
  }

  ## ---------------------------
  ## 7. Top markers heatmap
  ## ---------------------------
  top_markers <- unique(top_markers[!is.na(top_markers)])
  write.csv(
    top_markers,
    file      = file.path(out_dir, "Top_markers.csv"),
    row.names = FALSE,
    quote     = FALSE
  )
  message("\nTotal unique top markers: ", length(top_markers))

  if (length(top_markers) > 0) {
    message("Generating heatmap of top markers...")
    top_markers_valid <- intersect(top_markers, rownames(seurat[[assay]]))

    if (length(top_markers_valid) == 0) {
      message("No valid top markers present in assay '", assay, "' - skipping heatmap.")
    } else {
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

      avg_exp <- avg_exp_list[[assay]]

      paletteLength <- 50
      myColor  <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(
        seq(-2, 0, length.out = ceiling(paletteLength / 2) + 1),
        seq(2 / paletteLength, 2, length.out = floor(paletteLength / 2))
      )

      pdf(file.path(out_dir, "Top_markers_HM.pdf"), height = 25, width = 10)
      print(
        pheatmap(avg_exp,
                 cluster_rows  = TRUE,
                 scale         = "row",
                 color         = myColor,
                 breaks        = myBreaks,
                 border_color  = NA,
                 main          = "Scaled average expression of top markers",
                 fontsize      = 6)
      )
      dev.off()
    }
  } else {
    message("No top markers to plot")
  }

  ## ---------------------------
  ## 8. Timing and return
  ## ---------------------------
  end <- Sys.time()
  message("\n=== Analysis Complete ===")
  message("Start time: ", start)
  message("End time:   ", end)
  message("Duration:   ", round(difftime(end, start, units = "mins"), 2), " minutes")
  message("Comparison: positive log2FC = higher in ", conditions[1])

  invisible(list(
    all_results     = all_results,
    summary         = summary_results,
    top_markers     = top_markers,
    design_formula  = design_formula,
    reduced_formula = if (test_type == "LRT") reduced_formula else NULL,
    conditions      = conditions,
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
