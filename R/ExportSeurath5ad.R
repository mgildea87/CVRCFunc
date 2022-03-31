#' Export integrated seurat assay for RNA velocity and CellRank analysis
#' @param seurat Path to Seurat .rds file
#' @param sample_ident Identity for samples. Genereally 'orig.ident'
#' @param dir Directory to output h5ad file
#' @return h5ad file with integrated assay
#' @import Seurat SeuratDisk
#' @export


ExportSeurath5ad <- function(seurat, sample_ident, dir){

  seurat_obj <- readRDS(file = seurat)

  cells <- Cells(seurat_obj)
  cells <- strsplit(cells, split = "_")
  cells_vec <- vector()
  for(i in cells){
      cells_vec <- c(cells_vec, substr(dplyr::last(i), start = 1, stop = nchar(dplyr::last(i))-2))
  }
  cells_vec <- paste(seurat_obj@meta.data[[sample_ident]],":",cells_vec,"x", sep = "")
  #Rename cells in scvelo format
  seurat_obj <- RenameCells(seurat_obj, new.names = cells_vec)
  #correct for stupid bug in seurat
  for(i in 1:length(seurat_obj@assays$SCT@SCTModel.list)){
    slot(seurat_obj$SCT@SCTModel.list[[i]], 'median_umi') = median(seurat_obj$SCT@SCTModel.list[[i]]@cell.attributes$umi)
  }
  slot(seurat_obj$integrated@SCTModel.list[[1]], 'median_umi') = median(seurat_obj$integrated@SCTModel.list[[1]]@cell.attributes$umi)
  SaveH5Seurat(seurat_obj, filename = paste(dir,"Seurat", sep = ""))
  Convert(paste(dir,"Seurat.h5seurat",sep = ""), dest = "h5ad", overwrite = T, assay = 'integrated')
}