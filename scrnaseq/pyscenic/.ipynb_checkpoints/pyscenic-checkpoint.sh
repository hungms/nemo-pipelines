#!/bin/bash -l

#SBATCH --job-name=pyscenic
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

conda activate /camp/home/hungm/.conda/envs/pyscenic

cd /camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/anndata

python /camp/home/hungm/working/Matthew/project/anqi/20230719_SC22272/pipeline/16_pyscenic/pyscenic.py
