library(tidyverse)
library(Seurat)
library(SeuratData)
library(patchwork)
library(cowplot)
library(clustree)
library(slingshot)
library(tradeSeq)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e12)
setwd("/nemo/lab/caladod/working/Matthew/project/anqi/20230719_SC22272")

load("data/14_slingshot/14_slingshot_objects.Rdata")

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 30

set.seed(5)
cells_keep1 <- meta.data %>% filter(!is.na(C05_pseudotime)) %>% rownames(.)
cells_keep2 <- meta.data %>% filter(!is.na(C08_pseudotime)) %>% rownames(.)
cells_keep3 <- meta.data %>% filter(!is.na(C07_pseudotime)) %>% rownames(.)

counts <- as.matrix(counts.matrix[features, c(cells_keep1, cells_keep2, cells_keep3)])

png("figures/14_velocity_trajectory_png/14.7_slingshot/between_lineage_knots.png", width = 1500*5, height = 500*5, res = 78*5)
icMat <- evaluateK(counts = counts, pseudotime = pseudotime[cells_keep, c(6,3,1)], cellWeights = cell.weights[cells_keep, c(6,3,1)], k = 3:20, nGenes = 200, verbose = T, parallel=TRUE, BPPARAM = BPPARAM)
dev.off()
