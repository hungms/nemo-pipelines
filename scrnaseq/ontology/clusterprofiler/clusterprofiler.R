.libPaths("/nemo/lab/caladod/working/Matthew/.conda/envs/scrna/lib/R/library")
# load R libraries
library(tidyverse)
library(patchwork)
library(cowplot)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(readxl)
library(writexl)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e12)
setwd("/camp/home/hungm/scratch/hungm/RN23365/results")

# hallmark gene set
HM<-msigdbr(species ="Homo sapiens", category ="H")
HMgene<-HM %>% select(., gs_name, gene_symbol) %>% mutate(gs_name = gsub("HALLMARK_", "",gs_name))

# go gene set
GO<-msigdbr(species ="Homo sapiens", category ="C5")
GOgene<-GO %>% select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'GOBP')) %>% mutate(gs_name = gsub("GOBP_", "",gs_name))

gmt <- read.csv("/camp/home/hungm/scratch/hungm/RN23365/results/ssgsea/gene_sets.gmt", header = F, sep = "\t", row.names = 1)
gmt <- gmt[,c(3:ncol(gmt))]
DNA_damage_gene <- gmt %>% filter(rownames(.) == "GOBP_DNA_DAMAGE_RESPONSE") %>% pivot_longer(everything(), names_to = "gs_name", values_to = "gene_symbol")
DNA_damage_gene$gs_name <- rep("DNA_DAMAGE_RESPONSE", nrow(DNA_damage_gene))
GOgene <- rbind(GOgene, DNA_damage_gene)

# kegg gene set
KEGG<-msigdbr(species ="Homo sapiens", category ="C2")
KEGGgene<-KEGG %>% select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'KEGG')) %>% mutate(gs_name = gsub("KEGG_", "",gs_name))

# reactome gene set
REACTOME<-msigdbr(species ="Homo sapiens", category ="C2")
REACTOMEgene<-REACTOME %>% select(., gs_name, gene_symbol) %>% filter(stringr::str_detect(gs_name, 'REACTOME')) %>% mutate(gs_name = gsub("REACTOME_", "",gs_name))


ichop_chop_list <- list()
for(t in c("6h", "12h", "24h")){
    ichop_chop <- read.csv(paste0("C5_ichop_vs_chop_at_", t, ".csv"))
    ichop_chop_list[[length(ichop_chop_list) + 1]] <- ichop_chop}
names(ichop_chop_list) <- c("ichop_chop_6h", "ichop_chop_12h", "ichop_chop_24h")

for(x in c(1:3)){
    deglist <- ichop_chop_list[[x]] %>% arrange(log2FoldChange) %>% .$log2FoldChange
    names(deglist) <- ichop_chop_list[[x]] %>% arrange(log2FoldChange) %>% .$V3
    deglist <- na.omit(deglist)
    deglist <- deglist[which(deglist != 0)]
    deglist = sort(deglist, decreasing = TRUE)

    assign(paste0(names(ichop_chop_list)[x], "_HM"), GSEA(deglist, nPermSimple = 10000, TERM2GENE = HMgene, pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 500))
    assign(paste0(names(ichop_chop_list)[x], "_GO"), GSEA(deglist, nPermSimple = 10000, TERM2GENE = GOgene, pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 500))
    assign(paste0(names(ichop_chop_list)[x], "_KEGG"), GSEA(deglist, nPermSimple = 10000, TERM2GENE = KEGGgene, pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 500))
    assign(paste0(names(ichop_chop_list)[x], "_REACTOME"), GSEA(deglist, nPermSimple = 10000, TERM2GENE = REACTOMEgene, pvalueCutoff = 1.1, pAdjustMethod = "fdr", minGSSize = 1, maxGSSize = 500))
}

save(list = ls(pattern="ichop_chop_.*h_"), file = "/camp/home/hungm/scratch/hungm/RN23365/results/gsea/ichop_chop_GSEA.Rdata")

