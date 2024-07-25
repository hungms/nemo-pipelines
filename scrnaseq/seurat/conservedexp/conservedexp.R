source("/camp/home/hungm/working/Matthew/library/rlib/singlecell/seurat/seurat.R")
setwd("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/")

#rna <- qread("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
#rna@meta.data$celltypes <- ifelse(as.character(rna@meta.data$celltypes) == "MYC+CCDN3+", "MYC+CCND3+", as.character(rna@meta.data$celltypes))
#rna@meta.data$celltypes <- factor(rna@meta.data$celltypes, c("Activated", "GC-DZ", "GC-LZ", "Memory", "MYC+CCND3+"))
#rna@meta.data$leiden_0.6 <- factor(rna@meta.data$leiden_0.6, 0:12)
#rna@meta.data$group_labels <- gsub("GFPpos", "", rna@meta.data$group)
#rna@meta.data$group_labels <- gsub("B220pos", "PLACEBO", rna@meta.data$group_labels)
#metadata <- rna@meta.data %>% 
#    rownames_to_column("cell_barcode") %>%
#    group_by(celltypes, leiden_0.6) %>%
#    mutate(celltypes_leiden = as.factor(paste0(celltypes, "_", leiden_0.6))) %>%
#    ungroup() %>%
#    column_to_rownames("cell_barcode")
#rna@meta.data <- metadata[colnames(rna),]


#rna <- NormalizeData(rna)
#Idents(rna) <- "leiden_0.6"
#levels <- levels(rna@meta.data$leiden_0.6)

#markers <- list()
#for(x in 1:length(levels)){
#  markers[[x]] <- FindConservedMarkers(
#    rna,
#    ident.1 = levels[x],
#    grouping.var = "group_labels",
#    assay = "RNA",
#    slot = "data",
#    method = "MAST",
#    min.cells.group = 3)
  
#  markers[[x]] <- markers[[x]] %>%
#    rownames_to_column("gene") %>%
#    mutate(ident = levels[x])
#}

#qsave(markers, "output/conservedexp/20240717_integrated_conservedexp.qs")

markers <- qread("output/conservedexp/20240717_integrated_conservedexp.qs")
cols <- unique(unlist(lapply(markers, colnames)))

add_missing_cols <- function(df, all_cols) {
  missing_cols <- setdiff(all_cols, colnames(df))
  df[missing_cols] <- NA
  return(df)}
markers <- lapply(markers, add_missing_cols, all_cols = cols)

# Bind rows
markers <- bind_rows(markers)
write.csv(markers, "output/conservedexp/20240717_integrated_conservedexp.csv", quote = F, row.names = T)
