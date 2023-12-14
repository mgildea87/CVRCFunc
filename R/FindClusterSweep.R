#' Run common single-cell RNA-seq clustering algorithms as implemented in Seurat across a range of resolution values and compute common clustering metrics. This function assumes a KNN graph already exists in the specified assay. Run seurat's FindNeighbors before this function.
#' @param seurat A seurat object
#' @param assay seurat assay e.g. 'RNA'
#' @param resolutions Vector of clustering resolutions
#' @param algorithm Seurat FindClusters algorithm parameter. From Seurat: 'Algorithm for modularity optimization 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm. Leiden requires the leidenalg python.'
#' @param conda If applicable, path to conda environment. Required for Leiden algorithm = 4
#' @param pca_dim vector of dimensions to use for computing silheoutte scores. default = 0,30
#' @param reduction reduction to use for computing silhouette scores. default = 'pca'
#' @param plot_reduction reduction to plot silhouette scores. default = 'umap'
#' @return A seurat object with clustering. .pdf document with a series of clustering related plots
#' @import Seurat clustree ggplot2 bluster cluster pheatmap igraph patchwork
#' @importFrom reticulate use_miniconda
#' @importFrom ggbeeswarm geom_quasirandom
#' @export


FindClusterSweep <- function(seurat, assay = 'RNA', resolutions = c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6), algorithm = 1, conda, pca_dim = c(1:30), reduction = 'pca', plot_reduction = 'umap'){
  #Clustering
  DefaultAssay(seurat) <- assay
  if(algorithm == 4){
    use_miniconda(condaenv = conda)
    for(i in resolutions){
      print(paste0('clustering with resolution = ',i,'   ',Sys.time()))
      seurat <- FindClusters(seurat, resolution = i, algorithm = algorithm, method = "igraph")
    }
  }
  else{
    for(i in resolutions){
      print(paste0('clustering with resolution = ',i,'   ',Sys.time()))
      seurat <- FindClusters(seurat, resolution = i, algorithm = algorithm)
    }
  }
  print(paste0('Plotting...',Sys.time()))

  #Plotting
  #cluster numbers
  cluster_idents <- paste0(assay,'_snn_res.',resolutions)
  cluster_numbers <- vector()
  for(i in cluster_idents){
    cluster_numbers <- c(cluster_numbers, length(unique(seurat@meta.data[,which(colnames(seurat@meta.data) == i)])))
  }
  gg_frame_cluster_number <- data.frame(resolution = resolutions, clusters = as.numeric(cluster_numbers))
  pdf('FindClusterSweep_plots.pdf')
  print(ggplot(gg_frame_cluster_number, aes(y=clusters,x=resolution))+
    geom_line(color = 'steelblue')+
    ylab('Number of clusters')+
    xlab('resolution')+
    theme_bw())
  print(clustree::clustree(seurat, prefix = paste0(assay,"_snn_res.")))

  #silhouette scores
  dist.matrix <- dist(x = Embeddings(object = seurat[[reduction]])[, pca_dim])
  silhouette_scores <- data.frame(row.names = colnames(seurat))
  for(i in cluster_idents){
    RNA_sil <- silhouette(x = as.numeric(x = as.factor(x = seurat@meta.data[[i]])), dist = dist.matrix)
    RNA_sil <- as.data.frame(RNA_sil)
    RNA_sil$cluster <- as.character(RNA_sil$cluster)
    seurat[[paste0('silscore_',i)]] <- RNA_sil$sil_width
    colnames(RNA_sil) <- paste0(i,"_",colnames(RNA_sil))
    silhouette_scores <- cbind(silhouette_scores, RNA_sil)
  }
  silhouette_scores_scores <- silhouette_scores[,grep(colnames(silhouette_scores), pattern = 'sil_width')]
  colnames(silhouette_scores_scores) <- resolutions
  silhouette_scores_scores <- melt(silhouette_scores_scores)
  colnames(silhouette_scores_scores)[1] <- 'resolution'
  print(ggplot(silhouette_scores_scores, aes(y = value, x = resolution, color = resolution))+
    ggbeeswarm::geom_quasirandom()+
    geom_boxplot(alpha = 0.5, color = 'black', outlier.shape = NA)+
    theme_bw()+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
    ylab('Silhouette scores')+
    xlab('Resolution'))

  #per resolution plots
  plots_sil_emb <- list()
  plots_emb <- list()
  for(i in cluster_idents){
    RNA_sil <- silhouette_scores[,grep(colnames(silhouette_scores), pattern = i)]
    colnames(RNA_sil) <- c('cluster', 'neighbor', 'silhouette_score')
    RNA_sil$closest <- factor(ifelse(RNA_sil$silhouette_score > 0, RNA_sil$cluster, RNA_sil$neighbor))
    print(ggplot(RNA_sil, aes(y = silhouette_score, x = cluster, color = closest))+
      geom_quasirandom(alpha = 0.6, shape=16)+
      theme_classic()+
      ylab('Silhouette width')+
      ggtitle(i))

    gg_frame_silwidth <- data.frame(umap1 = seurat@reductions[[plot_reduction]]@cell.embeddings[,1], umap2 = seurat@reductions[[plot_reduction]]@cell.embeddings[,2], sil_width = seurat[[paste0('silscore_',i)]][,1])
    p_sil_emb <- ggplot(gg_frame_silwidth, aes(x = umap1, y = umap2, color = sil_width))+
      geom_point()+
      scale_color_gradient2(low="blue", mid="lightgrey", high="red", limits = c(-1,1))+
      theme_classic()+
      ggtitle(i,'\ncells colored by silhouette width')
    print(p_sil_emb)
    plots_sil_emb[[i]] <- ggplot(gg_frame_silwidth, aes(x = umap1, y = umap2, color = sil_width))+geom_point(size = .1)+
      scale_color_gradient2(low="blue", mid="lightgrey", high="red", limits = c(-1,1))+
      theme_classic()+theme(text = element_text(size=4))+ggtitle(i,'\ncells colored by silhouette width')

    mod <- pairwiseModularity(igraph::graph_from_adjacency_matrix(seurat@graphs[[paste0(assay,'_nn')]]), seurat[[i]][,1], as.ratio = TRUE)
    print(pheatmap(log10(mod+1), cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("white", "red"))(100), display_numbers = T, main = paste0('Pairwise cluster modularity\n',i)))
    cluster.gr <- igraph::graph_from_adjacency_matrix(log2(mod+1),mode="upper", weighted=TRUE, diag=FALSE)
    set.seed(11001010)
    print(plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5, layout=igraph::layout_with_lgl))

    Idents(seurat) <- i
    p_emb <- DimPlot(seurat, label = T, pt.size = .1, label.size = 2)+NoLegend()+ggtitle(i)
    print(p_emb)
    plots_emb[[i]] <- p_emb+theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(), text = element_text(size=4))

  }
  grid.arrange(grobs = plots_sil_emb, ncol = 2)
  print(wrap_plots(plots_emb, ncol = 3))
  dev.off()
  print(Sys.time())
  return(seurat)
}
