#' Differential expression testing between conditions via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param condition_ident Identity class for conditions to be tested
#' @param conditions A vector of the 2 conditions within condition_ident to be included in the DESeq2 model.
#' @param expfilt genes that have greater than 0 counts in greater than expfilt fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @return .csv files with marker genes per clus_ident. .pdf files with plots
#' @import Seurat pheatmap DESeq2 Matrix.utils reshape2 ggplot2 ggrepel stringr utils grDevices BiocGenerics
#' @importFrom BiocGenerics t
#' @importFrom magrittr set_colnames
#' @export

FindMarkersCondition <- function(seurat, clus_ident, sample_ident, condition_ident, conditions, expfilt = 0.5){
  start <- Sys.time()

  dir.create("FindMarkersCondition_outs", showWarnings = FALSE)

  Idents(seurat) <- clus_ident
  clusters <- unique(Idents(seurat))

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
  pb <- split.data.frame(pb, factor(splitf)) %>% lapply(function(u) set_colnames(t(u), str_extract(rownames(u), "(?<=_)[:alnum:]+")))

  # Print out the table of cells in each cluster-sample group
  pdf(paste('FindMarkersCondition_outs/cells_per_clus_HM.pdf'))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  for(cluster in clusters){

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
    print(DESeq2::plotPCA(vst, intgroup = "sample") +theme_classic() +geom_text_repel(aes(label = sample)))
    plotDispEsts(dds)
    plotMA(res)
    dev.off()
    #Save results table
    res <- as.data.frame(res_shrink)
    write.csv(file = paste('FindMarkersCondition_outs/cluster_',cluster,"_results.csv", sep = ''), res)
  }
  print(start)
  print(Sys.time())
}

