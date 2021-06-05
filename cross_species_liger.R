#load packages
library(liger, lib.loc = "/storage/chen/home/qingnanl/anaconda3/envs/bindSC/lib/R/library/")
library(Seurat)
library(dplyr)
library(Matrix)
library(cluster)
library(SeuratWrappers)
library(ape)

#BC
BC_species <- readRDS("~/BC_species1.rds")#this is a Seurat object with human, mouse and monkey BCs, the key was 'species', and there is another
#key named 'cell_identity', which is like "human_BC0" etc..


#feature selection
ifnb.list <- SplitObject(BC_species, split.by = "species")
ifnb.list.1 <- lapply(X = ifnb.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "dispersion", nfeatures = 2000)
  })
features = Reduce(union, list(VariableFeatures(ifnb.list.1[[1]]), VariableFeatures(ifnb.list.1[[2]]), VariableFeatures(ifnb.list.1[[3]])))


#run liger
BC_species <- ScaleData(BC_species, features = features, split.by = "species", do.center = FALSE)
BC_species <- RunOptimizeALS(BC_species, k = 40, lambda = 5, split.by = "species")#
BC_species <- RunQuantileNorm(BC_species, split.by = "species")
BC_species <- FindNeighbors(BC_species, reduction = "iNMF", dims = 1:40)
BC_species <- FindClusters(BC_species, resolution = 0.1)
BC_species <- RunUMAP(BC_species, dims = 1:ncol(BC_species[["iNMF"]]), reduction = "iNMF")


#plot to see the integration result
tiff("./liger_2/BC_species_40k.tiff", res = 300, width = 3000, height = 2000)
p<-DimPlot(BC_species, group.by = "species")
print(p)
dev.off()


#use 40 reduced dims to calculate distance
BC_dim<-as.data.frame(t(BC_species@reductions$iNMF@cell.embeddings[, 1:40]))
#average by cell_identity
BC_dim_avg <- sapply(split(colnames(BC_dim), as.character(BC_species@meta.data$cell_identity)),
                         function(cells) rowMeans(as.matrix(BC_dim[,cells])))
BC_dim_avg<-as.data.frame(t(BC_dim_avg))
dist <- dist(BC_dim_avg , diag=TRUE)
phy1 <- nj(dist)#tree file

#plot the tree
pdf("~/BC_species_1_tree_liger.pdf")
p<-plot(phy1, edge.width = 2)
print(p)
dev.off()

#now need to calculate a consensus tree, because we don't want the result to be affected by parameter selection that much
treeList<-c()


for (j in 1:100){
  df<-data.frame(cell_identity = BC_species@meta.data$cell_identity, index = 1:length(BC_species@meta.data$cell_identity))
  new_df <- df %>% group_by(cell_identity) %>% sample_n(100)
  sel<-new_df$index
  BC <- BC_dim[, sel]#select columns
  #print(head(BC))
  BC<-BC[sample(c(1:nrow(BC)), floor(0.8*nrow(BC)), replace = F), ]
  #print(head(BC))
  BC_int_avg <- sapply(split(colnames(BC), new_df$cell_identity), function(cells) rowMeans(as.matrix(BC[,cells])))
  dist <- dist(as.data.frame(t(BC_int_avg)), diag=TRUE)
  phy1 <- nj(dist)
  #nam <- paste("tree", j, sep = "")
  #assign(nam, phy1)
  treeList[j]<-c(phy1)#pass tree to vector
  print(j)
}

tree<-consensus(treeList, p = 1, check.labels = TRUE)
pdf("~/BC_species_consensus.pdf")
p<-plot(tree, edge.width = 2)
print(p)
dev.off()






