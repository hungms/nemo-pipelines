#!/bin/bash
#SBATCH --job-name=fastq_zip
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:0
#SBATCH --mem=100G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

####### edit here ########
export accession=GSE219098
#########################

mkdir -p /camp/home/hungm/scratch/hungm/${accession}/raw_fastq/tmp
cd /camp/home/hungm/scratch/hungm/${accession}/raw_fastq/tmp
mv *.gz ..
mv ../*.fastq .
conda activate jupyterlab
ls | xargs -n 1 -P 30 pigz
