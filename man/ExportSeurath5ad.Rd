% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportSeurath5ad.R
\name{ExportSeurath5ad}
\alias{ExportSeurath5ad}
\title{Export integrated seurat assay for RNA velocity, CellRank, Scanpy, and other python based analysis}
\usage{
ExportSeurath5ad(seurat, sample_ident, dir = "", assay)
}
\arguments{
\item{seurat}{Path to Seurat .rds file}

\item{sample_ident}{Identity for samples. Genereally 'orig.ident'}

\item{dir}{Directory to output h5ad file}

\item{assay}{which assay to save}
}
\value{
h5ad file with integrated assay
}
\description{
Export integrated seurat assay for RNA velocity, CellRank, Scanpy, and other python based analysis
}
