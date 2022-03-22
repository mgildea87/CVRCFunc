#' Export meta data files for RNA velocity and CellRank analysis
#' @param seurat Path to Seurat .rds file
#' @param clus_ident Identity for clusters. Normally 'seurat_clusters' but can be any identity. I generally use 'annotations' to store cluster annotations
#' @param sample_ident Identity for samples. Genereally 'orig.ident'
#' @param dir Directory to output .csv files
#' @return .csv files with cell by embeddings and clusters. .csv file with color palette for plotting in scVelo, CellRank, etc.
#' @import Seurat
#' @importFrom scales hue_pal
#' @export


ExportSeuratMeta <- function(seurat, clus_ident, sample_ident, dir){

  seurat_obj <- readRDS(file = seurat)

  cells <- Cells(seurat_obj)
  cells <- strsplit(cells, split = "_")
  cells_vec <- vector()
  for(i in cells){
      cells_vec <- c(cells_vec, substr(dplyr::last(i), start = 1, stop = nchar(dplyr::last(i))-2))
  }
  cells_vec <- paste(seurat_obj@meta.data[[sample_ident]],":",cells_vec,"x", sep = "")
  write.table(cells_vec, file = paste(dir,"cellID_obs.csv",sep = ""), row.names = FALSE, quote = T)

  emb_frame <- data.frame(row.names = cells_vec)
  for(i in 1:length(seurat_obj@reductions)){
    emb_frame <- cbind(emb_frame, Embeddings(seurat_obj@reductions[[i]])[,1:2])
  }
  col_name <- vector()
  for(i in names(CD4@reductions)){
    col_name <- c(col_name, paste(toupper(i),"_1", sep = ""), paste(toupper(i),"_2", sep = ""))
  }
  colnames(emb_frame) <- col_name
  write.csv(emb_frame, file = paste(dir,"cell_embeddings.csv",sep=""), row.names = T)

  clusters <- data.frame(clusters = seurat_obj@meta.data[[clus_ident]])
  row.names(clusters) <- cells_vec
  write.csv(clusters, file = paste(dir,clus_ident,".csv", sep = ""), row.names = T)

  colors <- hue_pal()(length(unique(seurat_obj@meta.data[[clus_ident]])))
  write.csv(colors, file = paste(dir,clus_ident,"_colors.csv", sep = ""), row.names = F)
}
