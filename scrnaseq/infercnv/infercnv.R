###
source("/camp/home/hungm/working/Matthew/library/rlib/singlecell/seurat/seurat.R")
library(infercnv)
library(qs)
setwd("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/output/infercnv")
gene_annot <- read.table("/flask/scratch/caladod/hungm/reference/alignment/infercnv/GRCm39_gene_annotations.txt", sep = "\t", row.names = 1, header = F)
colnames(gene_annot) <- NULL

SC23423 = qread("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
reference = qread("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/seurat/20240628_SC23423_B220posGFPneg.qs")
genes <- intersect(rownames(SC23423@assays$RNA$counts), rownames(reference@assays$RNA$counts))

SC23423[["RNA_filter"]] <- CreateAssay5Object(SC23423@assays$RNA$counts[genes,], min.cells = 5)
SC23423@meta.data$group <- gsub("GFPpos", "", SC23423@meta.data$group)
SC23423@meta.data$group	<- gsub("B220pos", "PLACEBO", SC23423@meta.data$group)
SC23423@meta.data$leiden_0.6 = paste0("C", sprintf("%02d", as.numeric(SC23423@meta.data$leiden_0.6)-1))
print(unique(SC23423@meta.data$leiden_0.6))

reference[["RNA_filter"]] <- CreateAssay5Object(reference@assays$RNA$counts[genes,], min.cells = 5)
reference@meta.data$group <- "B220+GFP-"
reference@meta.data$leiden_0.6 <- "B220+GFP-"

SC23423_total <- merge(SC23423, reference)
SC23423_total[["RNA_filter"]] <- JoinLayers(SC23423_total[["RNA_filter"]])
rna_counts = SC23423_total@assays$RNA_filter$counts
annotations = SC23423_total@meta.data %>%
		select("group")
colnames(annotations) <- NULL


rm(SC23423, reference)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=rna_counts,
                                    annotations_file=annotations,
                                    delim="\t",
                                    gene_order_file=gene_annot,
                                    ref_group_names=c("B220+GFP-"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/infercnv/20240705_group-level/",
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE,
			     HMM_type = c("i6"),
			     output_format="pdf",
			     num_threads = 30
)

qsave(infercnv_obj, file = "infercnv_obj.qs")
