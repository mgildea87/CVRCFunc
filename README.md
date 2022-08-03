# CVRCFunc
To install run `devtools::install_github('mgildea87/CVRCFunc')`
This is a simple package that will include R functions I think will be useful for others in the CVRC. The idea is to have them in a centralized location so we can all use them.

## FindMarkersBulk()
This function performs a similar task to Seurat's FindAllMarkers() but does so by pseudobulking by sample and then performing DESeq2. For each cluster, counts from cells are pseudobulked from the cluster of interest and from all other clusters. This generates 2 pseudobulks for each sample for each cluster. The cluster of interest pseudobulks are then compared to pseudobulks from all other clusters using DESeq2.

## FindMarkersCondition()
This function performs differential expression for each cluster (or other identity specified) between samples from 2 specified conditions. This is done via pseudobulking and DESeq2.

## ExportSeuratMeta()
This function exports seurat meta data for use in RNA velocity analysis. specifically, .csv files of reformatted sample/cell barcodes, cell by cluster mappings, cell by embedding coordinates, and color palette for plotting.

## ExportSeurath5ad()
Converts and saves a seurat object integrated assay as .h5ad. For use in RNA velocity analysis.

## StackedVlnPlot()
Creates a stacked violin plot from a Seurat object and a vector of gene names. Meant to emulate the Scanypy function.
