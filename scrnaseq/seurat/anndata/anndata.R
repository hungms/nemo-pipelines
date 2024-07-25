library(tidyverse)
library(Seurat)
library(SeuratData)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e13)

###### change everything contained here ######
setwd("/camp/home/hungm/working/Matthew/project/oscar/OA_ZND6629388")
load("seurat/ZND6629388_B.Rdata")
obj.name <- ls()
obj <- get(obj.name)
pca.name <- "pca"
reduction.name <- "umap"
##############################################

#obj.list <- SplitObject(obj, split.by = "sample.id")
#for(x in 1:length(obj.list)){#
#	obj.list[x] <- RenameCells(obj.list[x], add.cell.id = unique(obj@meta.data$sample.id))}
#obj <- merge(obj.list[[1]], obj.list[2:length(obj.list)])
#obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#obj[["SCT"]] <-	JoinLayers(obj[["SCT"]])

#names <- paste0(obj@meta.data$sample.id, "_", gsub("-1_.*$", "-1", colnames(obj)))
#obj <- RenameCells(obj, new.names = names)

#vf <- VariableFeatures(obj[["SCT"]])
#write.csv(vf, file = paste0("anndata/variable_features.csv"), quote = F)

#count
#obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
counts <- obj[["RNA.MM"]]$counts
write.csv(counts, file = paste0("anndata/rna_counts.csv"))
rm(counts)

#data
#data <- obj[["RNA.MM"]]$data
#write.csv(data, file = paste0("anndata/rna_data.csv"))
#rm(data)

#metadata
metadata <- obj@meta.data
write.csv(metadata, file = paste0("anndata/metadata.csv"))
rm(metadata)

#neighbour
#nn <- obj@graphs$SCT_nn
#write.csv(nn, file = paste0("anndata/nn.csv"))
#rm(nn)

#snn <- obj@graphs$SCT_snn
#write.csv(snn, file = paste0("anndata/snn.csv"))
#rm(snn)

#pca
pca <- obj@reductions[[paste0(pca.name)]][[]]
write.csv(pca, file = "anndata/pca.csv")
rm(pca)

#reductions
reduction <- obj@reductions[[paste0(reduction.name)]][[]]
write.csv(reduction, file = "anndata/reduction.csv")
rm(reduction)

rm(list = ls())
