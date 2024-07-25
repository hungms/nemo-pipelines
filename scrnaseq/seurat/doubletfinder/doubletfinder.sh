#!/bin/bash -l

#SBATCH --job-name=doubletfinder
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --time=3-00:00
#SBATCH --mem=200000
#SBATCH --partition=cpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

export DATA=/camp/home/hungm/working/Matthew/project/lingling/SC23355/
mkdir -p ${DATA}/logs/seurat
rm /camp/home/hungm/template_pips/seurat/doubletfinder/slurm-*
cp -r /camp/home/hungm/template_pips/seurat/doubletfinder/ ${DATA}logs/seurat/
exec > ${DATA}/logs/seurat/doubletfinder/doubletfinder.log 2>&1

conda activate /camp/home/hungm/.conda/envs/seurat5
/camp/home/hungm/.conda/envs/seurat5/bin/Rscript ${DATA}logs/seurat/doubletfinder/doubletfinder.R
