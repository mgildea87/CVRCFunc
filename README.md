# CVRCFunc
To install run `devtools::install_github('mgildea87/CVRCFunc')`
This is a simple package that will include R functions I think will be useful for others in the CVRC. The idea is to have them in a centralized location so we can all use them.

## FindMarkersBulk()
This function performs a similar task to Seurat's FindAllMarkers() but does so by pseudobulking by sample and then performing DESeq2. For each cluster, counts from cells are pseudobulked from the cluster of interest and from all other clusters. This generates 2 pseudobulks for each sample for each cluster. The cluster of interest pseudobulks are then compared to pseudobulks from all other clusters using DESeq2.
