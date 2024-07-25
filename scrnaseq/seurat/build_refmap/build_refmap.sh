#!/bin/bash -l

#SBATCH --job-name=refmap
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export logroot=/camp/home/hungm/working/Matthew/project/anqi/AX_GSE230705/logs/
mkdir -p ${logroot}
mkdir -p ${logroot}seurat/
mkdir -p ${logroot}seurat/build_refmap/

rm -r /camp/home/hungm/working/Matthew/template_pips/seurat/build_refmap/slurm-*
cp -r /camp/home/hungm/working/Matthew/template_pips/seurat/build_refmap/build_refmap* ${logroot}seurat/build_refmap/
exec > ${logroot}seurat/build_refmap/build_refmap.log 2>&1

conda activate /camp/home/hungm/.conda/envs/scrna
Rscript ${logroot}seurat/build_refmap/build_refmap.R
