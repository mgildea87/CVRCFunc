#' Differential expression testing between 2 specified groups of cells via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param group_1 Identity of cells within \code{clus_ident} to compare
#' @param group_2 Identity of cells within \code{clus_ident} to compare
#' @param expfilt_freq genes that have greater than \code{expfilt_counts} in greater than \code{expfilt_freq} fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @param expfilt_counts genes with less than \code{expfilt_counts} in \code{expfilt_freq * sample number} will be removed from DESeq2 model. 1 by default.
#' @param alpha FDR adjusted p-value threshold for significance in plotting. 0.05 by default.
#' @param out_dir Name of output directory
#' @return .csv files with marker genes per \code{clus_ident}. .pdf files with diagnostic plots
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom magrittr set_colnames
#' @export

FindMarkers <- function(seurat, clus_ident, group_1, group_2, sample_ident, expfilt_counts = 1, expfilt_freq = 0.5, out_dir = "FindMarkers", alpha = 0.05){
  start <- Sys.time()

  coef <- variable <- value <- NULL

  dir.create(out_dir, showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  seurat <- subset(seurat, idents = c(group_1, group_2))
  clusters <- unique(Idents(seurat))
  clusters <- sort(clusters)

  # Print out the table of cells in each cluster-sample group
  seurat@meta.data[,clus_ident] <- droplevels(factor(seurat@meta.data[,clus_ident]))
  pdf(paste(out_dir,'/',group_1,"_",group_2,"_",'cells_per_group_HM.pdf',sep = ''))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 8)
  dev.off()

  #Run wilcoxauc from presto for pct.in and pct.out (maybe save the stats later for comparison?)
  wilcox <- wilcoxauc(seurat, assay = 'data', seurat_assay = 'RNA', group_by = clus_ident)

  # Subset metadata to only include the cluster and sample IDs to aggregate across
  groups <- seurat@meta.data[, c(clus_ident, sample_ident)]

  # Aggregate across cluster-sample groups
  pb <- aggregate.Matrix(t(seurat@assays$RNA@counts), groupings = groups[,1:2], fun = "sum")

  # Not every cluster is present in all samples; create a vector that represents how to split samples
  splitf <- sapply(stringr::str_split(rownames(pb), pattern = "_",  n = 2), `[`, 2)

  cluster_counts <- as.data.frame(t(as.matrix(pb)))
  cluster_metadata <- data.frame(cluster = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 1),
                                 sample = sapply(stringr::str_split(colnames(cluster_counts), pattern = "_",  n = 2), `[`, 2))
  cluster_metadata$cluster <- factor(cluster_metadata$cluster, levels = c(group_1, group_2))

  dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ cluster)
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

  #Plots
  gg_counts <- cluster_counts[,sort(colnames(cluster_counts))]
  gg_counts <- melt(log10(gg_counts))
  pdf(paste(out_dir,"/",colnames(coef(dds))[2],"_diagnostic_plots.pdf", sep = ''))
  print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
  print(DESeq2::plotPCA(vst, intgroup = "cluster") +theme_classic())
  print(DESeq2::plotPCA(vst, intgroup = "sample") +theme_classic() +geom_text_repel(aes(label = sample), show.legend = FALSE))
  plotDispEsts(dds)
  print(ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig))+geom_point(alpha = 0.7)+scale_color_manual(values = c('grey40', 'blue'))+theme_classic())
  plotMA(res)
  dev.off()


  #Add pct.in and pct.out and save results table
  wilcox_groups <- unique(wilcox$group)
  wilcox_sub <- wilcox[which(wilcox$group == wilcox_groups[1]),]
  res_shrink$feature <- row.names(res_shrink)
  merged_res <- merge(x = res_shrink, y = wilcox_sub, by = 'feature', sort = F)
  res_shrink[,paste('pct_in',wilcox_groups[1],sep = "_")] <- merged_res$pct_in
  res_shrink[,paste('pct_in',wilcox_groups[2],sep = "_")] <- merged_res$pct_out
  res_shrink$feature <- NULL
  write.csv(file = paste(out_dir,"/",colnames(coef(dds))[2],"_results.csv", sep = ''), res_shrink)

  print(start)
  print(Sys.time())
}
