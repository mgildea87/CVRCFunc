#' Differential expression testing between conditions via pseudobulking and DESeq2
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity
#' @param sample_ident Sample identities. Identity class that indicates how to partition samples
#' @param condition_ident Identity class for conditions to be tested
#' @param conditions A vector of the 2 conditions within \code{condition_ident} to be included in the DESeq2 model.
#' @param expfilt_freq genes that have greater than \code{expfilt_counts} in greater than \code{expfilt_freq} fraction of cells will be kept for the DESeq2 model. 0.5 by default
#' @param expfilt_counts genes with less than \code{expfilt_counts} in \code{expfilt_freq * sample number} will be removed from DESeq2 model. 1 by default.
#' @param out_dir Name of output directory
#' @param alpha FDR adjusted p-value threshold for significance in plotting. 0.1 by default.
#' @return .csv files with marker genes per \code{clus_ident}. .pdf files with diagnostic plots
#' @import Seurat pheatmap DESeq2 reshape2 ggplot2 ggrepel stringr utils grDevices
#' @importFrom BiocGenerics t
#' @importFrom magrittr set_colnames
#' @export

FindMarkersCondition <- function(seurat, clus_ident, sample_ident, condition_ident, conditions, expfilt_counts = 1, expfilt_freq = 0.5, out_dir = "FindMarkersCondition_outs", alpha = 0.1){
  start <- Sys.time()

  coef <- variable <- value <- NULL

  dir.create(out_dir, showWarnings = FALSE)

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
  pdf(paste(out_dir,'/cells_per_clus_HM.pdf', sep = ''))
  pheatmap(table(seurat@meta.data[,clus_ident], seurat@meta.data[,sample_ident]), display_numbers = T, cluster_rows = F, cluster_cols = F, fontsize_number = 4)
  dev.off()

  sig_up <- vector()
  sig_down <- vector()
  clusters_summary <- vector()

  for(cluster in clusters){

    cond_1_cells <- which(seurat@meta.data[[clus_ident]] == cluster & seurat@meta.data[[condition_ident]] == conditions[1])
    cond_2_cells <- which(seurat@meta.data[[clus_ident]] == cluster & seurat@meta.data[[condition_ident]] == conditions[2])
    pct_in_cond_1 <- apply(seurat@assays$RNA@counts[,cond_1_cells,drop=F], MARGIN = 1, function (x) sum(x > 0) / length(x))
    pct_in_cond_2 <- apply(seurat@assays$RNA@counts[,cond_2_cells,drop=F], MARGIN = 1, function (x) sum(x > 0) / length(x))
    pct_in <- data.frame(pct_in_cond_1, pct_in_cond_2, feature = names(pct_in_cond_1))

    cluster_counts <- as.data.frame(as.matrix(pb[[as.character(cluster)]]))
    cluster_metadata <- samplecon_table[which(samplecon_table$sample %in% colnames(cluster_counts)),]
    cluster_metadata$condition <- factor(cluster_metadata$condition, levels = c(conditions[2], conditions[1]))
    cluster_counts <- cluster_counts[,cluster_metadata$sample]

    #skip loop iteration here if there are < 2 replicates in one or both conditions
    if(length(which(cluster_metadata$condition == levels(cluster_metadata$condition)[1])) < 2 | length(which(cluster_metadata$condition == levels(cluster_metadata$condition)[2])) < 2 ){
      print(paste('One or both conditions have fewer than 2 replicates for cluster:',cluster))
      next
    }

    dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ condition)
    # Filter data
    keep <- rowSums(counts(dds) >= expfilt_counts) >= expfilt_freq*nrow(cluster_metadata)
    dds <- dds[keep,]
    # DESeq2
    vst <- varianceStabilizingTransformation(dds)
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds, alpha = 0.1)

    clusters_summary <- c(clusters_summary, cluster)
    sig_up <- c(sig_up, length(which(res$padj < alpha & res$log2FoldChange > 0)))
    sig_down <- c(sig_down, length(which(res$padj < alpha & res$log2FoldChange < 0)))

    # Shrink the log2 fold changes to be more appropriate using the apeglm method
    res_shrink <- lfcShrink(dds, coef = colnames(coef(dds))[2], res=res, type = "apeglm")
    res_shrink <- as.data.frame(res_shrink)
    res_shrink$sig <- rep('Not significant', nrow(res_shrink))
    res_shrink$sig[which(res_shrink$padj < alpha)] <- 'Significant'

    #set up data for top gene plot
    d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
    exp_gene_logfc <- res[which.min(res$padj),]$log2FoldChange
    exp_gene_padj <- res[which.min(res$padj),]$padj

    #Plots
    gg_counts <- cluster_counts
    gg_counts <- melt(log10(gg_counts))
    pdf(paste(out_dir,'/cluster_',cluster,"_diagnostic_plots.pdf", sep = ''))
    print(ggplot(gg_counts, aes(x = variable, y = value, fill = variable)) + geom_boxplot() +  theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "none") + ylab('Log10(Counts)'))
    print(DESeq2::plotPCA(vst, intgroup = "condition") +theme_classic())
    print(DESeq2::plotPCA(vst, intgroup = "sample") +theme_classic() +geom_text_repel(aes(label = sample), show.legend = FALSE))
    plotDispEsts(dds)
    plotMA(res)
    print(ggplot(res_shrink, aes(x = log2FoldChange, y = -log10(pvalue), color = sig))+geom_point(alpha = 0.7)+scale_color_manual(values = c('grey40', 'blue'))+theme_classic())
    print(ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + ggtitle(paste('Gene:', row.names(res)[which.min(res$padj)],'\nLog2FC = ',exp_gene_logfc,'\npadj = ',exp_gene_padj, sep = '')))
    dev.off()

    #Save results table
    res_shrink$feature <- row.names(res_shrink)
    merged_res <- merge(x = res_shrink, y = pct_in, by = 'feature', sort = F)
    res_shrink[[paste('pct_in_',conditions[1], sep = '')]] <- merged_res$pct_in_cond_1
    res_shrink[[paste('pct_in_',conditions[2], sep = '')]] <- merged_res$pct_in_cond_2
    res_shrink$feature <- NULL

    write.csv(file = paste(out_dir,'/cluster_',cluster,"_results.csv", sep = ''), res_shrink)
  }
  summary_results <- data.frame(cluster = clusters_summary, sig_up = sig_up, sig_down = sig_down)
  write.csv(file = paste(out_dir,'/Summary_results.csv', sep = ''), summary_results, quote = F, row.names = F)
  print(start)
  print(paste('Order of comparison in DESeq2 model.   ', colnames(coef(dds))[2]))
  print(Sys.time())
}

