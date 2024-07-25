#!/bin/bash -l

#SBATCH --job-name=cibersortx
#SBATCH --partition=ncpu
#SBATCH --ntasks=1
#SBATCH --time=7-00:00
#SBATCH --mem=250000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### change data path ###
export accession=oscar/OA_SC23423
export DATA=
export LOGS=
export sc_counts= #full path
export bulk_counts= #full path

########################
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}/
source /camp/home/hungm/pipelines/piplog.sh

singularity exec \
-B $DATA:/src/data/ \
-B $LOGS:/src/outdir \
$sif /src/CIBERSORTxFractions \
--username matthew.hung@crick.ac.uk \
--token 0d87b21cb59616468e2b83abde26f771 \
--single_cell TRUE  \
--refsample $sc_counts \
--mixture $bulk_counts \
--verbose TRUE
