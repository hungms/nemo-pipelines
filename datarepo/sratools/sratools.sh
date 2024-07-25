#!/bin/bash
#SBATCH --job-name=sratools
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:0
#SBATCH --mem=200G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

############## edit the following ##################
export accession=anqi/AX_GSE219098
export PRJ=/camp/home/hungm/scratch/hungm/${accession}
#####################################################

#conda activate /camp/home/hungm/.conda/envs/sratools
source /camp/home/hungm/pipelines/piplog.sh

which fasterq-dump
mkdir -p ${PRJ}/raw_fastq/
cd ${PRJ}/raw_fastq/
for i in $(cat ${PRJ}/input/SRR_Acc_List.txt);
	do ~/.conda/envs/sratools/bin/fasterq-dump $i -e 32 --include-technical --split-files;
	for j in i ; 
		do gzip ${i}*.fastq ; 
	done; 
done
