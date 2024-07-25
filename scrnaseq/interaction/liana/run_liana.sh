#!/bin/bash -l

#SBATCH --job-name=liana
#SBATCH --partition=cpu
#SBATCH --ntasks=21
#SBATCH --time=3-00:00
#SBATCH --mem=250000
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

########################
### change data path ###
export DATA=/camp/home/hungm/working/Matthew/project/lingling/LZ_SC23355/
########################

mkdir -p ${DATA}/logs/interaction
rm /camp/home/hungm/template_pips/interaction/liana/slurm-*
cp -r /camp/home/hungm/template_pips/interaction/liana/ ${DATA}logs/interaction/
exec > ${DATA}/logs/interaction/liana/${SLURM_JOB_ID}_run_liana.log 2>&1

conda activate /camp/home/hungm/.conda/envs/liana
python ${DATA}/logs/interaction/liana/run_liana.py
