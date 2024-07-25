#!/bin/bash -l

#SBATCH --job-name=pearson
#SBATCH --partition=ncpu
#SBATCH --ntasks=4
#SBATCH --time=1-00:00
#SBATCH --mem=20G

conda activate /camp/home/hungm/.conda/envs/seurat5
/camp/home/hungm/.conda/envs/seurat5/bin/Rscript /camp/home/hungm/pipelines/tmp/tmp.R
