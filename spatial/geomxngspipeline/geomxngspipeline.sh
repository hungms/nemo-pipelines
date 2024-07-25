#!/bin/sh
#SBATCH --job-name=geomxngspipeline
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=%j.out
#SBATCH --error=%j.err

GEOMXNGS_HOME=/scratch/prj/cb_lymphnode/geomxpipeline/2.3.3.10
DATASET_HOME=/scratch/prj/cb_lymphnode/Matthew/AC_S947
RESULTS_HOME=/scratch/prj/cb_lymphnode/Matthew/AC_S947

$GEOMXNGS_HOME/geomxngspipeline -in $DATASET_HOME/fastq --out=$RESULTS_HOME/dcc --ini=$DATASET_HOME/ini/231221_AC_GeomX_run2_20240108T1613_GNP_config.ini --save-interim-files=true --threads=3
