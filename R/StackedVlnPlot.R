#' Generate a stacked violin plot from a Seurat object and vector of genes.
#' @param seurat A Seurat object
#' @param clus_ident Identity for clusters. 'seurat_clusters' by default.
#' @param features vector of gene names to plot
#' @param assay assay to plot from. default = 'RNA'
#' @param slot vslot to plot from. default = 'data'
#' @param main Title for plot
#' @import Seurat patchwork ggplot2 purrr
#' @export


StackedVlnPlot <- function(seurat, clus_ident = 'seurat_clusters', features, assay = 'RNA', slot = 'data', main = ''){
  Idents(seurat) <- clus_ident
  DefaultAssay(seurat) <- assay

  modify_vlnplot<- function(obj,
                            feature,
                            pt.size = 0,
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {
    p<- VlnPlot(obj, features = feature, pt.size = pt.size, slot = slot, ... )  +
      xlab("") + ylab(feature) + ggtitle("") +
      theme(legend.position = "none",
            plot.title= element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = rel(1), angle = 0),
            axis.text.y = element_text(size = rel(1)),
            plot.margin = plot.margin )
    return(p)
  }

  ## extract the max value y-axis value
  extract_max<- function(p){
    ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
    return(ceiling(ymax))
  }


  Plot<- function(obj, features,
                            pt.size = 0,
                            plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                            ...) {

    plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

    plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
      theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1), axis.ticks.x = element_line())

    # change the y-axis tick to max value
    ymaxs<- purrr::map_dbl(plot_list, extract_max)
    plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                              scale_y_continuous(breaks = c(y)) +
                              expand_limits(y = y))

    p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) + patchwork::plot_annotation(title = main, theme = theme(plot.title = element_text(size = 18)))
    return(p)
  }
  return(Plot(obj = seurat, features = features))
}

