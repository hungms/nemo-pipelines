library(tidyverse)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(clustree)
options(Seurat.object.assay.version = "v4")
options(future.globals.maxSize = 1e12)

### change here ###
setwd("/nemo/lab/caladod/working/Matthew/project/anqi/AX_GSE230705")
dataset <- "GSE230705"
split <- "sample.id"
###################

counts <- read.csv("anndata/rna_counts.csv", row.names = 1)
colnames(counts) <- gsub("X", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))

v4 <- CreateSeuratObject(counts = counts)
metadata <- read.csv("anndata/metadata.csv", row.names = 1)
v4@meta.data <- metadata
all_features <- rownames(v4)
v4.list <- SplitObject(v4, split.by = split)

for(x in 1:length(v4.list)){
  v4.list[[x]] <- SCTransform(v4.list[[x]], vst.flavor = "v2", variable.features.n = 3000, assay = "RNA", return.only.var.genes = FALSE)
  v4.list[[x]] <- RunPCA(v4.list[[x]], npcs = 50, verbose = FALSE, assay = "SCT")}

features <- read.csv("anndata/variable_features.csv", sep = ",")$x
v4.list <- PrepSCTIntegration(object.list = v4.list, anchor.features = features)
v4.anchors <- FindIntegrationAnchors(object.list = v4.list, normalization.method = "SCT", anchor.features = features, reduction = "rpca", k.filter = 80)
v4.integrated <- IntegrateData(anchorset = v4.anchors, normalization.method = "SCT", features.to.integrate = all_features, k.weight = 50)
save(v4.integrated, file = paste0("seurat/", dataset, "_v4RefMap.Rdata"))

data <- as.matrix(v4.integrated[["integrated"]]@data)
write.csv(data, file = "anndata/integrated_data.csv", row.names = T)
scale.data <- as.matrix(v4.integrated[["integrated"]]@scale.data)
write.csv(scale.data, file = "anndata/integrated_scale.csv", row.names = T)

