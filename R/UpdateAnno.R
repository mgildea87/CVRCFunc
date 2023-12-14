#' Update cell annotations.
#' @param seurat Seurat object
#' @param convert_file absolute path to .csv file containing 2 columns. 1st column should have the header identical to \code{old_anno_ident} and contain each annotation to update as a row. The 2nd column should have a header identical to \code{new_anno_ident} and contain the updated annotations matched to \code{old_anno_ident}
#' @param old_anno_ident Identity of annotations to be changed
#' @param new_anno_ident Identity for new annotations
#' @return Seurat object with updated cell annotations specified in \code{convert_file}
#' @import Seurat
#' @export


UpdateAnno <- function(seurat, old_anno_ident, new_anno_ident, convert_file){
  convert_file <- read.csv(convert_file)
  new_annotations <- as.vector(seurat[[old_anno_ident]][,1])
  for(i in 1:nrow(convert_file)){
    new_annotations[which(new_annotations == convert_file[i,1])] <- convert_file[i,2]
  }
  seurat[[new_anno_ident]] <- new_annotations
  return(seurat)
}
