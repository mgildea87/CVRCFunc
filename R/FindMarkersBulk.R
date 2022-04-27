#' Find all markers via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param expfilt genes that have greater than 0 counts in greater than expfilt fraction of cells will be kept for the DESeq2 model
#' @param n_top_genes number of top genes per cluster to save and make a heatmap with
#' @param pct.in Filter threshold for top marker genes per cluster. For a given gene and cluster, if the fraction of cells with counts exceeds pct.in it is removed from top_markers.
#' @return .csv files with marker genes per clus_ident. .pdf files with plots
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom presto wilcoxauc
#' @export

FindMarkersBulk <- function(seurat, clus_ident, sample_ident, expfilt = 0.5, n_top_genes = 50, pct.in = 25){
  start <- Sys.time()

  coef <- variable <- value <- NULL

  dir.create("FindMarkersBulk_outs", showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  # Print out the table of cells in each cluster-sample group
  pdf(paste('FindMarkersBulk_outs/cells_per_clus_HM.pdf'))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  #Run wilcoxauc from presto for pct.in and pct.out (maybe save the stats later for comparison?)
  wilcox <- wilcoxauc(seurat, assay = 'data', seurat_assay = 'RNA', group_by = clus_ident)

  top_markers <- vector()

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
    vst <- varianceStabilizingTransformation(dds)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, alpha = 0.05)

    # Shrink the log2 fold changes to be more appropriate using the apeglm method
    res_shrink <- lfcShrink(dds, coef = colnames(coef(dds))[2], res=res, type = "apeglm")

    #Plots
    gg_counts <- cluster_counts[,sort(colnames(cluster_counts))]
    gg_counts <- melt(log10(gg_counts))
    pdf(paste('FindMarkersBulk_outs/cluster_',cluster,"_diagnostic_plots.pdf", sep = ''))
    print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
    print(DESeq2::plotPCA(vst, intgroup = "iscluster") +theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample_ident") +theme_classic() +geom_text_repel(aes(label = sample_ident), show.legend = FALSE))
    plotDispEsts(dds)
    plotMA(res)
    dev.off()


    #Save results table
    res_shrink <- as.data.frame(res_shrink)

    #Add pct.in and pct.out
    wilcox_sub <- wilcox[which(wilcox$group == cluster),]
    res_shrink$feature <- row.names(res_shrink)
    merged_res <- merge(x = res_shrink, y = wilcox_sub, by = 'feature', sort = F)
    res_shrink$pct_in <- merged_res$pct_in
    res_shrink$pct_out <- merged_res$pct_out
    res_shrink$feature <- NULL
    write.csv(file = paste('FindMarkersBulk_outs/cluster_',cluster,"_results.csv", sep = ''), res_shrink)

    res_sig <- res_shrink[which(res_shrink$padj < 0.05 & res_shrink$pct_in > pct.in),]
    top_markers <- c(top_markers,row.names(res_sig[order(res_sig$log2FoldChange),])[1:n_top_genes])
  }
  write.csv(top_markers, file = 'FindMarkersBulk_outs/Top_markers.csv', row.names = F, quote = F)
  pdf(file = 'FindMarkersBulk_outs/Top_markers_HM.pdf')
  print(DoHeatmap(subset(seurat, downsample = 1000), features = top_markers, assay = 'RNA', slot = 'scale.data', raster = F)+
          scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill", na.value = "white") +
          theme(text = element_text(size = 1)))
  dev.off()
  print(start)
  print(Sys.time())
}
