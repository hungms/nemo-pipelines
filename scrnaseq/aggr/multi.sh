#!/bin/bash
#SBATCH --job-name=cellranger-multi
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:0
#SBATCH --mem=250G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

export accession=anqi/AX_GSE219098
export sample=
mkdir -p /camp/home/hungm/scratch/hungm/${accession}/cellranger
exec > /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${sample}/multi.log 2>&1
rm -r /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/slurm*

cd /camp/home/hungm/scratch/hungm/${accession}/cellranger
module load CellRanger/7.1.0
cellranger multi --id=${sample} \
		 --csv=/camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${sample}/config.csv \
                 --localmem=200 \
                 --localcores=30

