#!/bin/bash -l

#SBATCH --job-name=singleseqgset
#SBATCH --partition=cpu
#SBATCH --ntasks=4
#SBATCH --time=3-00:00
#SBATCH --mem=200000

export logroot=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/

mkdir -p ${logroot}functional
mkdir -p ${logroot}functional/biological_processes

rm -r /camp/home/hungm/working/Matthew/template_pips/functional/biological_processes/singleseqgset/slurm-*
cp -r /camp/home/hungm/working/Matthew/template_pips/functional/biological_processes/singleseqgset/ ${logroot}functional/biological_processes
exec > ${logroot}functional/biological_processes/singleseqgset/singleseqgset.log 2>&1

conda activate /camp/home/hungm/.conda/envs/seurat5
Rscript ${logroot}functional/biological_processes/singleseqgset/singleseqgset.R
