#!/bin/bash -l

#SBATCH --job-name=dandelion_tcr
#SBATCH --partition=ncpu
#SBATCH --ntasks=15
#SBATCH --time=3-00:00
#SBATCH --mem=150000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

#----------------------------------------------
### specify samples and files - please change
export accession=/lingling/LZ_SC23355/

#----------------------------------------------
### environment setups - do not need to change
ml Singularity/3.11.3
export PRJ=/flask/scratch/caladod/hungm/$accession/
source /camp/home/hungm/pipelines/piplog.sh
export DATA=$PRJ/dandelion/tcr/
cd $DATA # change path

### keep logs
export LOGS=/flask/scratch/caladod/hungm/$accession/logs/scvdjseq/dandelion_updated/tcr/
export sif=/camp/home/hungm/pipelines/containers/sc-dandelion_latest.sif # set singularity container

#----------------------------------------------
### run singularity - please check and change
singularity run -B $DATA,$LOGS $sif dandelion-preprocess \
	--org mouse \
	--chain TR \
	--db imgt \
	--file_prefix all \
	--meta $LOGS/tcr_meta.csv \
	--keep_trailing_hyphen_number \
	--skip_tigger


# optional arguments
        #--org #human #mouse
        #--chain #IG #TR
        #--db #imgt #ogrdb (pick ogrdb if specify mouse strain)
        #--file_prefix #all #filter
        #--meta $SCRIPT/$meta
        #--keep_trailing_hyphen_number
        #--skip_tigger
        #--skip_reassign_dj
        #--skip_correct_c
