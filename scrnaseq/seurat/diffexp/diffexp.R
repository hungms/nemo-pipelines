source("/camp/home/hungm/working/Matthew/library/rlib/singlecell/seurat/seurat.R")
setwd("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/")

SC23423 <- qread("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
DefaultAssay(SC23423) <- "RNA"
SC23423	<- NormalizeData(SC23423)

SC23423$group_label <- gsub("GFPpos", "", SC23423$group)
SC23423$group_label <- gsub("B220pos", "PLACEBO", SC23423$group_label)


SC23423 <- subset(SC23423, subset = group == "TUMOUR")
Idents(SC23423) <- "leiden_0.6"

markers <- FindAllMarkers(SC23423,
                          test.use = "MAST",
                          only.pos = FALSE,
                          min.pct = 0.1,
                          logfc.threshold = -Inf,
                          return.thresh = 1,
                          assay = "RNA")
write.csv(markers, file = "output/diffexp/20240708_SC23423_integrated_TUMOUR_leiden_0.6.csv")
