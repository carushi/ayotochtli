library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackVlnPlotGG <- function(data, features, cluster, header, count_data=FALSE) {
    matrix = data.frame(as.matrix(data))
    rownames(matrix) = 1:dim(matrix)[1]
    colnames(matrix) = features
    print(head(matrix))
    matrix <- matrix[,colSums(matrix) > 0]
    cluster <- cluster[rowSums(matrix) > 0]
    matrix <- matrix[rowSums(matrix) > 0,]
    if (dim(matrix)[1] == 0) return()
    sfeatures = colnames(matrix)
    matrix$Cell = rownames(matrix)
    matrix$Idents = cluster
    print(head(matrix))
    df <- reshape2::melt(matrix, id.vars = c("Cell", "Idents"), measure.vars = sfeatures,
                        variable.name = "Feat", value.name = "Expr")    
    df[,'Idents'] <- factor(df[,'Idents'])
    print(head(df))
    if (count_data) scale = 'area'
    else scale = 'width'
    g <- ggplot(df, aes(x=Idents, y=Expr, fill=Feat))+geom_violin(scale=scale, adjust=1.3)+facet_grid(rows=vars(Feat), scales = "free", switch = "y") + theme_cowplot(font_size = 12)
    g <- g+theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
        ggtitle("identities on x-axis") + xlab("Identity") + ylab("Expression Level")
    g <- g+scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
                           c(rep(x = "", times = length(x)-2), round(x[length(x) - 1], 3), "")) 
    pdf(paste0('violin_', header, '.pdf'), width=6, height=12)
    plot(g)
    dev.off()
    print(dim(df[df$Expr > 0,]))
    print(table(df[(df$Expr > 0) & (df$Feat == 'CD4'),c('Idents')]))
    g <- ggplot(df[df$Expr > 0,], aes(x=Idents, y=Expr, fill=Feat))+geom_violin(scale='width', adjust=1.3)+facet_grid(rows = vars(Feat), scales = "free_y", switch = "y") +  theme_cowplot(font_size = 12) 
    g <- g+theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.background = element_rect(fill = NA, color = "black"),
              strip.background = element_blank(),
              strip.text = element_text(face = "bold"),
              strip.text.y.left = element_text(angle = 0)) +
        ggtitle("identities on x-axis") + xlab("Identity") + ylab("Expression Level")
    g <- g+scale_y_continuous(expand = c(0, 0), position="right", labels = function(x) { if (length(x) == 1) return(round(x, 3));
                           c(rep(x = "", times = max(2, length(x))-2), round(x[max(1, length(x) - 1)], 3), "")}) 
    pdf(paste0('violin_', header, '_wo0.pdf'), width=6, height=12)
    plot(g)
    dev.off()
}