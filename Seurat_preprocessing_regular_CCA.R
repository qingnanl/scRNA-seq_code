library(Seurat)
library(Matrix)
library(dplyr)
head<-"" #in the dir to file, there is a place that folder is named by sample name/chip name etc. head is the dir before it.
tail<-"/outs/filtered_feature_bc_matrix/"
setwd("~/")#set to working directory
sample_list<-list.files(getwd())
i<-1
data.list<-list()
for (sample in sample_list){
  address<-paste0(head, sample, tail)
  data <- Read10X(data.dir = address)
  data <- CreateSeuratObject(counts = data,  min.cells = 0.005*length(colnames(data)), min.features = 1000) #cutoff can be changed
  data@meta.data$sample<-sample #create metadata ahead
  data@meta.data$condition<-gsub("\\_.*","",sample) #this is specific. It has the WT/KO information in the sample name
  data.list[i]<-data
  i<-i+1
  }
names(data.list) <- sample_list
#below is to use original preprocessing pipeline: normalize, find variable features, then integrate.
data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
#CCA
data.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:20)
data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:20)
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
# For non-SCT, Scale is required
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 20, verbose = FALSE)
saveRDS(data.combined, "data.combined.rds")
# t-SNE/UMAP and Clustering
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:15, umap.method='umap-learn',metric = 'correlation')

data.combined <- RunTSNE(data.combined, reduction = "pca", dims = 1:15)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:15)
data.combined <- FindClusters(data.combined, resolution = 0.1)
library(cowplot)
p1 <- DimPlot(data.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE)
pdf("overview.pdf", width = 10, height = 5)
plot_grid(p1, p2)
dev.off()
saveRDS(data.combined, "data.combined.rds")
