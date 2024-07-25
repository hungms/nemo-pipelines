library(tidyverse)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(DoubletFinder)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e13)
setwd("/camp/home/hungm/working/Matthew/project/matthew/GSE230705/")

load("seurat/GSE230705_PCA1.Rdata")
GSE230705.list <- SplitObject(GSE230705, split.by = "sample.id")
for(x in 1:length(GSE230705.list)){
    
      #take object
      obj <- GSE230705.list[[x]]

      #find pk
      pk <- paramSweep(obj, PCs = 1:20, sct = T)
      z <- summarizeSweep(pk, GT = FALSE)
      bcmvn <- as.data.frame(find.pK(z))
      bcmvn_y <- bcmvn %>% top_n(1, BCmetric)

      #print bcmvn_y
      bcmvn_df <- bcmvn[,c("pK","BCmetric")]
      names(bcmvn_df)[which(names(bcmvn_df) == "BCmetric")] <- names(GSE230705.list)[x]

      #add cluster annotations
      annotations <- obj@meta.data$seurat_clusters
      homotypic.prop <- modelHomotypic(annotations)

      #calculate doublet estimate - set doubletrate as 0.08
      nExp_poi <- round(0.08*nrow(obj@meta.data))
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

      #doublet finder
      obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = as.numeric(paste((bcmvn_y[,"pK"][[1]]))), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
      obj@meta.data %>% select(starts_with("pANN_0.25")) %>% colnames() -> pann
      obj <- doubletFinder(obj, PCs = 1:20, pN = 0.25, pK = bcmvn_y[,"pK"][[1]], nExp = nExp_poi.adj, reuse.pANN = pann , sct = T)

      #save doublet column
      obj@meta.data$doublet <- obj@meta.data[[ncol(obj@meta.data)]]

      #clean column
      obj@meta.data[,c(which(colnames(obj@meta.data) %>% str_detect("pANN|DF.class")))] <- NULL

      #go back to GSE230705.list
      GSE230705.list[[x]]@meta.data <- obj@meta.data
      message(paste0("Complete doubletfinder for batch ", x))
      }
GSE230705.list <- merge(GSE230705.list[[1]], GSE230705.list[2:length(GSE230705.list)])
GSE230705@meta.data <- GSE230705.list@meta.data
save(GSE230705, file = "seurat/GSE230705_Doublet.Rdata")


###################
#remember to change .sh file
###################
