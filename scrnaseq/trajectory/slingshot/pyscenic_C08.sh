#!/bin/bash -l

#SBATCH --job-name=scenic_ss_C08
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

conda activate /camp/home/hungm/.conda/envs/slingshot
Rscript /camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/pipeline/14_velocity_trajectory/slingshot/pyscenic_C08.R
