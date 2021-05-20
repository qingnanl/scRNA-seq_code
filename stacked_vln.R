#Seurat 4 already has the 'stack' option
#this can be used for Seurat 3
#modified from online codes
library(ggplot2)
modify_vlnplot<- function(obj,
feature,
pt.size = 0,
plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
...) {
p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
xlab("") + ylab(feature) + ggtitle("") +
theme(legend.position = "none",
axis.text.x = element_blank(),
#axis.text.y = element_blank(),
axis.ticks.x = element_blank(),
axis.ticks.y = element_line(),
axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
plot.margin = plot.margin )
return(p)
}

## main function
StackedVlnPlot<- function(obj, features,
pt.size = 0.0001,
plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
...) {

plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
theme(axis.text.x=element_text(), axis.ticks.x = element_line())

p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
return(p)
}
#to use
StackedVlnPlot(data, pt.size=0, features=marker_set)#data is a seurat object; marker_set is a vector of genes
