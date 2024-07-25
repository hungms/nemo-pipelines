source("/camp/home/hungm/working/Matthew/library/rlib/singlecell/seurat/seurat.R")
setwd("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/")

rna <- qread("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
rna <- NormalizeData(rna)
rna[["RNA"]] <- CreateAssay5Object(counts = rna[["RNA"]]$counts, data = rna[["RNA"]]$data, min.feature = 0, min.cell = 5)
rna$run_leiden <- paste0(rna$run_id, "_C", rna$leiden_0.6)
run_leiden <- AverageExpression(rna, group.by = "run_leiden")
run_leiden <- run_leiden$RNA %>%
    as.data.frame(.) %>%
    rownames_to_column("gene") %>%
    pivot_longer(!c("gene"), names_to = "run_leiden", values_to = "exprs") %>%
    mutate(run_leiden = gsub("-", "_", run_leiden))
head(run_leiden)
abundance <- rna@meta.data %>%
    group_by(run_id, leiden_0.6, run_leiden) %>%
    summarize(count = n()) %>%
    group_by(run_id) %>%
    mutate(pct = count*100/sum(count)) %>%
    ungroup() %>%
    dplyr::select(run_leiden, leiden_0.6, pct) %>%
    merge(., run_leiden, by = "run_leiden", all = T)
head(abundance)

df <- abundance
markers <- read.csv("output/diffexp/20240708_SC23423_integrated_TUMOUR_leiden_0.6.csv", row.names = 1)
markers$diff <- markers$pct.1 - markers$pct.2
markers.total <- markers %>% 
    filter(!str_detect(gene, "^Mt|^Rp[sl]|^H2-|^Ig[hkl][aecmgvdj]|^Gm[0-9]")) %>%
    filter(avg_log2FC > 0.5 & diff > 0 & p_val_adj < 0.05) %>%
    filter(pct.1 > 0.1) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>% 
    ungroup()

output <- list()
for(x in c(0:12)){
    coeff_list <- c()
    pvalue_list <- c()
    gene_list <- markers.total %>%
        filter(cluster == x) %>%
        .$gene

    if(length(gene_list) == 0){
        next}

    for(g in gene_list){

        coeff_random <- c()
        for(y in c(0:12)[which(!c(0:12) %in% x)]){

            # a vector of population pearson correlation between selected gene and other clusters
            other_cluster <- df %>% filter(gene == g & leiden_0.6 == y) %>% arrange(run_leiden) %>% .$pct
            other_geneexp <- df %>% filter(gene == g & leiden_0.6 == y) %>% arrange(run_leiden) %>% .$exprs
            coeff <- cor(other_cluster, other_geneexp)
            if(is.na(coeff)){coeff <- 0}
            coeff_random <- c(coeff_random, coeff)}
            
        # pearson correlation between selected gene and selected cluster
        selected_cluster <- df %>% filter(gene == g & leiden_0.6 == x) %>% arrange(run_leiden) %>% .$pct
        selected_geneexp <- df %>% filter(gene == g & leiden_0.6 == x) %>% arrange(run_leiden) %>% .$exprs
        coeff <- cor(selected_cluster, selected_geneexp)
        coeff_list <- c(coeff_list, coeff)
            
        # one-sample z-test - to see if observed correlation is significantly different to population correlation
        if(is.na(coeff)){coeff <- 0}
        pvalue <- BSDA::z.test(coeff_random, mu = coeff, sigma.x = sd(coeff_random))$p.value
        pvalue_list <- c(pvalue_list, pvalue)
    }

    cluster <- rep(x, length(gene_list))
    output[[length(output) + 1]] <- data.frame(cluster, gene_list, coeff_list, pvalue_list)
}

output <- bind_rows(output)
write.csv(output, file = "output/cluster_signature/gene_abundance/20240717_SC23423_TUMOUR_gene_abundance.csv")
