#' Generate violin plots from a Seurat object and vector of genes.
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. 'seurat_clusters' by default.
#' @param features vector of gene/feature names to plot
#' @param assay assay to plot from. default = 'RNA'
#' @param seurat seurat object
#' @param condition_ident name of metadata column to split plots by
#' @param idents_to_plot vector of a subset of specific idents within clus_ident to plot
#' @param ncol number of columns if plotting multiple genes
#' @import Seurat ggplot2 dplyr reshape2 ggbeeswarm
#' @export

Better_vlnplot <- function(seurat, clus_ident = 'seurat_clusters', condition_ident = vector(), features, assay = 'RNA', idents_to_plot = vector(), ncol = 1){

  if(length(idents_to_plot) == 0){
    idents <- unique(pull(seurat@meta.data, clus_ident))
  } else
    idents <- idents_to_plot

  if(length(condition_ident) == 0){
    meta_sub <- seurat@meta.data[which(pull(seurat@meta.data, clus_ident) %in% idents),c(clus_ident), drop=F]
    exp <- seurat@assays[[assay]]@data[which(row.names(seurat@assays[[assay]]@data) %in% features),row.names(meta_sub)]
    meta_sub <- cbind(meta_sub, t(exp))
    row.names(meta_sub) <- NULL
    meta_sub <- melt(meta_sub)
    meta_sub_wo0 <- meta_sub
    meta_sub_wo0$value[which(meta_sub_wo0$value == 0)] <- NA
    p <- ggplot(meta_sub, aes(x = .data[[clus_ident]], y = value, fill = .data[[clus_ident]]))+
      geom_quasirandom(data = meta_sub, aes(x = .data[[clus_ident]], y = value, fill = .data[[clus_ident]], alpha = I(ifelse(value == 0, 0, .5))), dodge.width = 1, size = 1, bandwidth = 20, width = .25, shape = 21, color = 'black')+
      geom_violin(alpha = 0.85, width = .5, color = 'grey20', scale = 'width')+
      facet_wrap(~variable, ncol = ncol)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 1), axis.title.x = element_blank(), panel.grid = element_blank())+
      ylab('log10(Normalized expression)')+
      theme(strip.background = element_blank(),strip.text = element_text(hjust = 0, face="bold"))
  }else{
    meta_sub <- seurat@meta.data[which(pull(seurat@meta.data, clus_ident) %in% idents),c(clus_ident, condition_ident)]
    exp <- seurat@assays[[assay]]@data[which(row.names(seurat@assays[[assay]]@data) %in% features),row.names(meta_sub)]
    meta_sub <- cbind(meta_sub, t(exp))
    row.names(meta_sub) <- NULL
    meta_sub <- melt(meta_sub)
    meta_sub_wo0 <- meta_sub
    meta_sub_wo0$value[which(meta_sub_wo0$value == 0)] <- NA
    p <- ggplot(meta_sub, aes(x = .data[[clus_ident]], y = value, fill = .data[[condition_ident]]))+
      geom_quasirandom(data = meta_sub, aes(x = .data[[clus_ident]], y = value, fill = .data[[condition_ident]], alpha = I(ifelse(value == 0, 0, .5))), dodge.width = .75, size = 1, bandwidth = 20, width = .15, shape = 21, color = 'black')+
      geom_violin(width = .75, alpha = 0.85, color = 'grey20', scale = 'width')+
      facet_wrap(~variable, ncol = ncol)+
      theme_bw()+
      theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 1), axis.title.x = element_blank(), panel.grid = element_blank())+
      ylab('log10(Normalized expression)')+
      theme(strip.background = element_blank(),strip.text = element_text(hjust = 0, face="bold"))
  }
  return(p)
}
