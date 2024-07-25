#!/bin/bash
#SBATCH --job-name=edgeturbo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:0
#SBATCH --mem=50G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

########## edit here ############
export accession=CRA004574
export PRJ=/camp/home/hungm/scratch/hungm/${accession}
#################################

source /camp/home/hungm/pipelines/piplog.sh
wget -c -r -np -k -L -p -P ${PRJ}  ftp://download.big.ac.cn/gsa2/${accession}

