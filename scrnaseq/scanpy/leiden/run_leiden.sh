#!/bin/bash
#SBATCH --job-name=leiden
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:0
#SBATCH --mem=200G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### edit the following  ###
export PRJ=/camp/home/hungm/working/Matthew/project/anqi/AX_GSE230705/
###########################

source /camp/home/hungm/pipelines/piplog.sh
conda activate /camp/home/hungm/.conda/envs/scanpy
python ${PRJ}/logs/scrnaseq/leiden/run_leiden.py
