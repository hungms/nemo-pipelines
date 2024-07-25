source("/camp/home/hungm/working/Matthew/library/rlib/seurat/scpa.R")
setwd("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/")

SC23423 <- qread("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
DefaultAssay(SC23423) <- "RNA"
SC23423	<- NormalizeData(SC23423)

for(x in seq_along(mm_gmt)){
    onetoall <- FindAllPathways(SC23423, pathways = mm_gmt[[x]], group = "leiden_0.6", parallel = T, core = 30)
    assign(paste0(names(mm_gmt)[x], "_onetoall"), onetoall)}
SC23423_leiden_0.6_onetoall <- mget(ls(pattern = "_onetoall"))
qsave(SC23423_leiden_0.6_onetoall, file = "/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/20240621_scpa/SC23423_leiden_0.6_comparisons.qs")
rm(SC23423_leiden_0.6_onetoall)

for(x in seq_along(mm_gmt)){
    onetoall <- FindAllPathways(SC23423, pathways = mm_gmt[[x]], group = "group", parallel = T, core = 30)
    assign(paste0(names(mm_gmt)[x], "_onetoall"), onetoall)}
SC23423_group_onetoall <- mget(ls(pattern = "_onetoall"))
qsave(SC23423_group_onetoall, file = "/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/20240621_scpa/SC23423_group_comparisons_RCHOPcontrol.qs")
