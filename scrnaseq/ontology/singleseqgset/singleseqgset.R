.libPaths("/nemo/lab/caladod/working/Matthew/.conda/envs/seurat5/lib/R/library")
library(tidyverse)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(msigdbr)
library(clusterProfiler)
library(singleseqgset)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e13)
setwd("/nemo/lab/caladod/working/Matthew/project/anqi/BMPC_timestamped")

load("finalised/seurat_clean.Rdata")

DefaultAssay(SC22272.PC) <- "RNA"
rna.data <- as(SC22272.PC[["RNA"]]$data, "sparseMatrix")
logfc <- logFC(cluster.ids=SC22272.PC@meta.data$Cluster.ID, expr.mat=rna.data)

hm.mouse <- msigdbr(species="Mus musculus",category="H")
hm.names <- unique(hm.mouse$gs_name)
hm.sets <- vector("list",length=length(hm.names))
names(hm.sets) <- hm.names
for (i in names(hm.sets)) {
    hm.sets[[i]] <- pull(hm.mouse[hm.mouse$gs_name==i,"gene_symbol"])}
hm.gsea <- wmw_gsea(expr.mat=rna.data, cluster.cells=logfc[[1]], log.fc.cluster=logfc[[2]], gene.sets=hm.sets)

go.mouse <- msigdbr(species="Mus musculus",category="C5")
go.names <- unique(go.mouse$gs_name)
go.sets <- vector("list",length=length(go.names))
names(go.sets) <- go.names
for (i in names(go.sets)) {
    go.sets[[i]] <- pull(go.mouse[go.mouse$gs_name==i,"gene_symbol"])}
go.sets <- go.sets[(which(str_detect(names(go.sets), "GO")))]
go.gsea <- wmw_gsea(expr.mat=rna.data, cluster.cells=logfc[[1]], log.fc.cluster=logfc[[2]], gene.sets=go.sets)


kegg.mouse <- msigdbr(species="Mus musculus",category="C2")
kegg.names <- unique(kegg.mouse$gs_name)
kegg.sets <- vector("list",length=length(kegg.names))
names(kegg.sets) <- kegg.names
for (i in names(kegg.sets)) {
    kegg.sets[[i]] <- pull(kegg.mouse[kegg.mouse$gs_name==i,"gene_symbol"])}
kegg.sets <- kegg.sets[(which(str_detect(names(kegg.sets), "KEGG")))]
kegg.gsea <- wmw_gsea(expr.mat=rna.data, cluster.cells=logfc[[1]], log.fc.cluster=logfc[[2]], gene.sets=kegg.sets)

reactome.mouse <- msigdbr(species="Mus musculus",category="C2")
reactome.names <- unique(reactome.mouse$gs_name)
reactome.sets <- vector("list",length=length(reactome.names))
names(reactome.sets) <- reactome.names
for (i in names(reactome.sets)) {
    reactome.sets[[i]] <- pull(reactome.mouse[reactome.mouse$gs_name==i,"gene_symbol"])}
head(names(reactome.sets)[(which(str_detect(names(reactome.sets), "REACTOME_")))])
reactome.gsea <- wmw_gsea(expr.mat=rna.data, cluster.cells=logfc[[1]], log.fc.cluster=logfc[[2]], gene.sets=reactome.sets)

save(list = ls(pattern = ".gsea"), file = "data/functional/biological_processes/singleseqgset/singleseqgset.Rdata")
