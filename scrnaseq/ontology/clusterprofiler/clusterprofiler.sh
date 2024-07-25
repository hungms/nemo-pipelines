#!/bin/bash -l

#SBATCH --job-name=clusterprofiler
#SBATCH --partition=ncpu
#SBATCH --ntasks=4
#SBATCH --time=3-00:00
#SBATCH --mem=200000

## remember to change logroot
export logroot=/camp/home/hungm/working/Matthew/project/oscar/OA_RN23365/logs/

mkdir -p ${logroot}functional
mkdir -p ${logroot}functional/ontology

rm -r /camp/home/hungm/working/Matthew/template_pips/functional/ontology/clusterprofiler/slurm-*
cp -r /camp/home/hungm/working/Matthew/template_pips/functional/ontology/clusterprofiler/ ${logroot}functional/ontology
exec > ${logroot}functional/ontology/clusterprofiler/clusterprofiler.log 2>&1

conda activate /camp/home/hungm/.conda/envs/scrna
Rscript ${logroot}functional/ontology/clusterprofiler/clusterprofiler.R
