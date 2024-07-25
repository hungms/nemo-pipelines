#!/bin/bash -l

#SBATCH --job-name=conservedexp
#SBATCH --partition=ncpu
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=250000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### change data path ###
export accession=oscar/OA_SC23423
########################
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}/
source /camp/home/hungm/pipelines/piplog.sh

conda activate /camp/home/hungm/.conda/envs/seurat5
/camp/home/hungm/.conda/envs/seurat5/bin/Rscript ${PRJ}logs/scrnaseq/seurat/conservedexp/${SLURM_JOB_ID}/conservedexp.R
