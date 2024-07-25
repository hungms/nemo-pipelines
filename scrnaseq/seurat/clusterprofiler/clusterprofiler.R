source("/camp/home/hungm/working/Matthew/library/rlib/rnaseq/clusterprofiler.R")
deg <- read.csv("/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/diffexp/20240627_SC23423_integrated_RCHOPvsPLACEBO.csv",  row.names = 1)

FindPathways(deg, "/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/clusterprofiler/20240627_SC23423_integrated_RCHOPvsPLACEBO_clusterprofiler.qs")
#FindAllPathways(deg, "/camp/home/hungm/working/Matthew/project/oscar/OA_SC23423/output/clusterprofiler/SC23423_integrated_leiden_0.6_clusterprofiler.qs")

