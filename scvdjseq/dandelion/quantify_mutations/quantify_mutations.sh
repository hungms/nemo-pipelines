#!/bin/bash -l

#SBATCH --job-name=quantify_mutations
#SBATCH --partition=ncpu
#SBATCH --ntasks=32
#SBATCH --time=7-00:00
#SBATCH --mem=1000000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

export accession=/anqi/AX_GSE230705
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}
source /camp/home/hungm/pipelines/piplog.sh

conda activate /camp/home/hungm/.conda/envs/scdandelion
python ${PRJ}/logs/scvdjseq/dandelion/quantify_mutations/quantify_mutations_template.py
