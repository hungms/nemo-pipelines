#!/bin/bash -l

#SBATCH --job-name=rmagic
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --time=3-00:00
#SBATCH --mem=200000

conda activate /camp/home/hungm/.conda/envs/seurat5
Rscript /camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/pipeline/rmagic.R
