#!/bin/bash -l

#SBATCH --job-name=infercnv
#SBATCH --partition=ncpu
#SBATCH --ntasks=30
#SBATCH --time=7-00:00
#SBATCH --mem=500000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### change data path ###
export accession=oscar/OA_SC23423
########################
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}/
source /camp/home/hungm/pipelines/piplog.sh

conda activate /camp/home/hungm/.conda/envs/seurat5
/camp/home/hungm/.conda/envs/seurat5/bin/Rscript ${PRJ}logs/scrnaseq/infercnv/${SLURM_JOB_ID}/infercnv.R
