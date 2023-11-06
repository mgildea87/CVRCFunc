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
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom presto wilcoxauc
#' @export

FindMarkersBulk <- function(seurat, clus_ident, sample_ident, expfilt_counts = 1, expfilt_freq = 0.5, n_top_genes = 50, pct.in = 25, out_dir = "FindMarkersBulk_outs", alpha = 0.1, assay = 'RNA'){
  start <- Sys.time()

  coef <- variable <- value <- NULL

  dir.create(out_dir, showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  # Print out the table of cells in each cluster-sample group
  pdf(paste(out_dir,'/cells_per_clus_HM.pdf',sep = ''))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  #Run wilcoxauc from presto for pct.in and pct.out (maybe save the stats later for comparison?)
  wilcox <- wilcoxauc(seurat, assay = 'data', seurat_assay = assay, group_by = clus_ident)

  top_markers <- vector()

  for(cluster in clusters){
    # Subset metadata to only include the cluster and sample IDs to aggregate across
    groups <- seurat@meta.data[, c(clus_ident, sample_ident)]
    groups$iscluster <- as.vector(groups[[clus_ident]])
    groups$iscluster[which(groups$iscluster != cluster)] <- "other"

    # Aggregate across cluster-sample groups
    pb <- aggregate.Matrix(t(seurat@assays[[assay]]@counts), groupings = groups[,2:3], fun = "sum")

    # Not every cluster is present in all samples; create a vector that represents how to split samples
    splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 2)
    splitf[which(splitf != cluster)] <- "other"

    cluster_counts <- as.data.frame(t(as.matrix(pb)))
    cluster_metadata <- data.frame(sample_ident = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 1),
                                   iscluster = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 2))
    cluster_metadata$iscluster <- factor(cluster_metadata$iscluster, levels = c('other', as.vector(cluster)))

    dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ iscluster)
    # Filter data
    keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq*nrow(cluster_metadata)
    dds <- dds[keep,]
    # DESeq2
    vst <- varianceStabilizingTransformation(dds)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, alpha = 0.05)

    # Shrink the log2 fold changes to be more appropriate using the apeglm method
    res_shrink <- lfcShrink(dds, coef = colnames(coef(dds))[2], res=res, type = "apeglm")
    res_shrink <- as.data.frame(res_shrink)
    res_shrink$sig <- rep('Not significant', nrow(res_shrink))
    res_shrink$sig[which(res_shrink$padj < alpha)] <- 'Significant'

    #set up data for top gene plot
    d <- plotCounts(dds, gene=which.min(res$padj), intgroup="iscluster", returnData=TRUE)
    exp_gene_logfc <- res[which.min(res$padj),]$log2FoldChange
    exp_gene_padj <- res[which.min(res$padj),]$padj

    #Plots
    gg_counts <- cluster_counts[,sort(colnames(cluster_counts))]
    gg_counts <- melt(log10(gg_counts))
    pdf(paste(out_dir,'/cluster_',cluster,"_diagnostic_plots.pdf", sep = ''))
    print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
    print(DESeq2::plotPCA(vst, intgroup = "iscluster") +theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample_ident") +theme_classic() +geom_text_repel(aes(label = sample_ident), show.legend = FALSE))
    plotDispEsts(dds)
    plotMA(res)
    print(ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig))+geom_point(alpha = 0.7)+scale_color_manual(values = c('grey40', 'blue'))+theme_classic())
    print(ggplot(d, aes(x=iscluster, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + ggtitle(paste('Gene:', row.names(res)[which.min(res$padj)],'\nLog2FC = ',exp_gene_logfc,'\npadj = ',exp_gene_padj, sep = '')))
    dev.off()


    #Save results table
    #Add pct.in and pct.out
    wilcox_sub <- wilcox[which(wilcox$group == cluster),]
    res_shrink$feature <- row.names(res_shrink)
    merged_res <- merge(x = res_shrink, y = wilcox_sub, by = 'feature', sort = F)
    res_shrink$pct_in <- merged_res$pct_in
    res_shrink$pct_out <- merged_res$pct_out
    res_shrink$feature <- NULL
    write.csv(file = paste(out_dir,'/cluster_',cluster,"_results.csv", sep = ''), res_shrink)

    res_sig <- res_shrink[which(res_shrink$padj < 0.05 & res_shrink$pct_in > pct.in),]
    top_markers <- c(top_markers,row.names(res_sig[order(res_sig$log2FoldChange, decreasing = T),])[1:n_top_genes])
  }
  write.csv(top_markers, file = paste(out_dir,'/Top_markers.csv', sep = ''), row.names = F, quote = F)
  avg_exp <- AverageExpression(seurat, assays = 'RNA', slot = 'data', group.by = clus_ident, features = top_markers)
  avg_exp <- avg_exp[[assay]]
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), seq(2/paletteLength, 2, length.out=floor(paletteLength/2)))
  pdf(file = paste(out_dir,'/Top_markers_HM.pdf',sep = ''), height = 25)
  print(pheatmap(avg_exp, cluster_rows = F, scale = 'row', color = myColor, breaks = myBreaks, border_color = NA, main = 'scaled average expression', fontsize = 6))
  dev.off()
  print(start)
  print(Sys.time())
}
