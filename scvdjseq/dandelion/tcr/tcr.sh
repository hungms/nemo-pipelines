#!/bin/bash -l

#SBATCH --job-name=dandelion_tcr
#SBATCH --partition=ncpu
#SBATCH --ntasks=15
#SBATCH --time=3-00:00
#SBATCH --mem=150000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk


### remember to change mouse/human
export accession=lingling/LZ_SC23355

### make sure dandelion directory is created
ml Singularity/3.11.3
export PRJ=/flask/scratch/caladod/hungm/${accession}
source /camp/home/hungm/pipelines/piplog.sh
export DATA=${PRJ}/dandelion/tcr/
mkdir -p ${DATA}

export LOGS=/flask/scratch/caladod/hungm/${accession}/logs/scvdjseq/dandelion/tcr/
mkdir -p ${LOGS}

### run singularity pipeline
export sif=/camp/home/hungm/working/Matthew/library/pythonlib/dandelion/sc-dandelion_latest.sif
cd ${DATA}
singularity run -B ${DATA},${LOGS} $sif python ${LOGS}/tcr.py \
	--org mouse \
	--file_prefix all \
	--chain TR \
	--meta ${LOGS}/meta.csv \
	--keep_trailing_hyphen_number \
	--skip_reassign_alleles



##   --org human \
##   --file_prefix filtered \
##   --chain TR \
##   --skip_reassign_alleles \
