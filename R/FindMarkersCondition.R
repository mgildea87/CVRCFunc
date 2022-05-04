#' Differential expression testing between conditions via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param condition_ident Identity class for conditions to be tested
#' @param conditions A vector of the 2 conditions within condition_ident to be included in the DESeq2 model.
#' @param expfilt genes that have greater than 0 counts in greater than expfilt fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @return .csv files with marker genes per clus_ident. .pdf files with plots
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom magrittr set_colnames
#' @export

FindMarkersCondition <- function(seurat, clus_ident, sample_ident, condition_ident, conditions, expfilt = 0.5){
  start <- Sys.time()

  coef <- variable <- value <- NULL

  dir.create("FindMarkersCondition_outs", showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  #Sample condition table
  samplecon_table <- unique(data.frame(sample = seurat@meta.data[[sample_ident]], condition = seurat@meta.data[[condition_ident]]))
  #Filter samplecon_table for conditions being tested
  samplecon_table <- samplecon_table[which(samplecon_table$condition %in% conditions),]

  #extract count matrix
  groups <- seurat@meta.data[, c(clus_ident, sample_ident)]
  pb <- aggregate.Matrix(t(seurat@assays$RNA@counts), groupings = groups, fun = "sum")
  split_sample <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 2)
  pb <- pb[which(split_sample %in% samplecon_table$sample),]
  splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 1)
  pb <- split.data.frame(pb, factor(splitf)) %>% lapply(function(u) set_colnames(t(u), str_extract(rownames(u), "(?<=_).+")))

  # Print out the table of cells in each cluster-sample group
  pdf(paste('FindMarkersCondition_outs/cells_per_clus_HM.pdf'))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  for(cluster in clusters){

    cond_1_cells <- which(seurat@meta.data[[clus_ident]] == cluster & seurat@meta.data[[condition_ident]] == conditions[1])
    cond_2_cells <- which(seurat@meta.data[[clus_ident]] == cluster & seurat@meta.data[[condition_ident]] == conditions[2])
    pct_in_cond_1 <- apply(seurat@assays$RNA@counts[,cond_1_cells,drop=F], MARGIN = 1, function (x) sum(x > 0) / length(x))
    pct_in_cond_2 <- apply(seurat@assays$RNA@counts[,cond_2_cells,drop=F], MARGIN = 1, function (x) sum(x > 0) / length(x))
    pct_in <- data.frame(pct_in_cond_1, pct_in_cond_2, feature = names(pct_in_cond_1))

    cluster_counts <- as.data.frame(as.matrix(pb[[as.character(cluster)]]))
    cluster_metadata <- samplecon_table[which(samplecon_table$sample %in% colnames(cluster_counts)),]
    cluster_counts <- cluster_counts[,cluster_metadata$sample]

    dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ condition)
    # Filter data
    keep <- rowSums(counts(dds) >= 1) >= expfilt*nrow(cluster_metadata)
    dds <- dds[keep,]
    # DESeq2
    vst <- varianceStabilizingTransformation(dds)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, alpha = 0.1)

    # Shrink the log2 fold changes to be more appropriate using the apeglm method
    res_shrink <- lfcShrink(dds, coef = colnames(coef(dds))[2], res=res, type = "apeglm")

    #Plots
    gg_counts <- cluster_counts
    gg_counts <- melt(log10(gg_counts))
    pdf(paste('FindMarkersCondition_outs/cluster_',cluster,"_diagnostic_plots.pdf", sep = ''))
    print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
    print(DESeq2::plotPCA(vst, intgroup = "condition") +theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample") +theme_classic() +geom_text_repel(aes(label = sample), show.legend = FALSE))
    plotDispEsts(dds)
    plotMA(res)
    dev.off()
    #Save results table
    res_shrink <- as.data.frame(res_shrink)

    res_shrink$feature <- row.names(res_shrink)
    merged_res <- merge(x = res_shrink, y = pct_in, by = 'feature', sort = F)
    res_shrink[[paste('pct_in_',conditions[1], sep = '')]] <- merged_res$pct_in_cond_1
    res_shrink[[paste('pct_in_',conditions[2], sep = '')]] <- merged_res$pct_in_cond_2
    res_shrink$feature <- NULL

    write.csv(file = paste('FindMarkersCondition_outs/cluster_',cluster,"_results.csv", sep = ''), res_shrink)
  }
  print(start)
  print(paste('Order of comparison in DESeq2 model.   ', colnames(coef(dds))[2]))
  print(Sys.time())
}

