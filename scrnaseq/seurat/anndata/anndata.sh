#!/bin/bash -l

#SBATCH --job-name=anndata
#SBATCH --partition=ncpu
#SBATCH --ntasks=4
#SBATCH --time=4-00:00
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

########################
### change data path ###
export accession=oscar/OA_ZND6629388

########################

export PRJ=/camp/home/hungm/working/Matthew/project/${accession}/
source /camp/home/hungm/pipelines/piplog.sh


conda activate /camp/home/hungm/.conda/envs/seurat5
/camp/home/hungm/.conda/envs/seurat5/bin/Rscript ${PRJ}/logs/scrnaseq/seurat/anndata/anndata.R

python ${PRJ}/logs/scrnaseq/seurat/anndata/anndata_save.py
