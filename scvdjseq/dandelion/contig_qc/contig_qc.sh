#!/bin/bash -l

#SBATCH --job-name=contigs_qc
#SBATCH --partition=ncpu
#SBATCH --ntasks=4
#SBATCH --time=3-00:00
#SBATCH --mem=50G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### set variables : accession, metadata, prj path, data path
export accession=/oscar/OA_SC23423/
export meta=/nemo/lab/caladod/working/Matthew/project/$accession/mudata/merged/SC23423_merged_rchop_placebo_tumour_integrated_metadata.csv
export PRJ=/flask/scratch/caladod/hungm/$accession/
export DATA=$PRJ/dandelion/20240627_bcr/
export sif=/camp/home/hungm/pipelines/scvdjseq/dandelion_updated/sc-dandelion_latest.sif #set singularity container

### setup to keep logs, specify scripts and singularity container - do not need to change
cp $meta /camp/home/hungm/pipelines/scvdjseq/dandelion_updated/contig_qc/metadata.csv
source /camp/home/hungm/pipelines/piplog.sh
ml Singularity/3.11.3
export LOGS=$PRJ/logs/scvdjseq/dandelion_updated/contig_qc/$now # make log file for contig_qc
cd $DATA
mkdir -p $DATA/contig_qc/ # make directory for contig_qc

### run singularity pipeline - please review
singularity run -B $DATA,$LOGS $sif python $LOGS/contig_qc.py \
	--meta $LOGS/metadata.csv \
	--chain IG \
	--file_prefix all \
	--keep_10x_c_call \
	--size 4

rm $LOGS/contig_qc.py
