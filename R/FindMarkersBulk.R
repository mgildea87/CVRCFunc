#' Fina all markers via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param expfilt genes that have greater than 0 counts in greater than expfilt fraction of cells will be kept for the DESeq2 model
#' @return .csv files with marker genes per clus_ident. .pdf files with plots
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices BiocGenerics
#' @importFrom BiocGenerics t
#' @export

FindMarkersBulk <- function(seurat, clus_ident, sample_ident, expfilt){
  start <- Sys.time()

  dir.create("FindMarkersBulk_outs", showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))

  # Print out the table of cells in each cluster-sample group
  pdf(paste('FindMarkersBulk_outs/cells_per_clus_HM.pdf'))
  pheatmap(table(seurat$seurat_clusters, seurat$new.ident), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  for(cluster in clusters){
    # Subset metadata to only include the cluster and sample IDs to aggregate across
    groups <- seurat@meta.data[, c(clus_ident, sample_ident)]
    groups$iscluster <- as.vector(groups[[clus_ident]])
    groups$iscluster[which(groups$iscluster != cluster)] <- "other"

    # Aggregate across cluster-sample groups
    pb <- aggregate.Matrix(t(seurat@assays$RNA@counts), groupings = groups[,2:3], fun = "sum")

    # Not every cluster is present in all samples; create a vector that represents how to split samples
    splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 2)
    splitf[which(splitf != cluster)] <- "other"

    cluster_counts <- as.data.frame(t(as.matrix(pb)))
    cluster_metadata <- data.frame(sample_ident = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 1),
                                   iscluster = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 2))

    dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ iscluster)
    # Filter data
    keep <- rowSums(counts(dds) >= 1) >= expfilt*nrow(cluster_metadata)
    dds <- dds[keep,]
    # DESeq2
    vst <- vst(dds, blind=TRUE)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, name = paste("iscluster_other_vs_",cluster, sep = ''), alpha = 0.1)

    # Shrink the log2 fold changes to be more appropriate using the apeglm method
    res <- lfcShrink(dds, coef = paste("iscluster_other_vs_",cluster, sep = ''), res=res, type = "apeglm")

    #Plots
    gg_counts <- cluster_counts[,sort(colnames(cluster_counts))]
    gg_counts <- melt(log10(gg_counts))
    pdf(paste('FindMarkersBulk_outs/cluster_',cluster,"_diagnostic_plots.pdf", sep = ''))
    print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
    print(DESeq2::plotPCA(vst, intgroup = "iscluster") +theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample_ident") +theme_classic() +geom_text_repel(aes(label = sample_ident)))
    plotDispEsts(dds)
    plotMA(res)
    dev.off()
    #Save results table
    res <- as.data.frame(res)
    write.csv(file = paste('FindMarkersBulk_outs/cluster_',cluster,"_results.csv", sep = ''), res)
  }
  print(start)
  print(Sys.time())
}
