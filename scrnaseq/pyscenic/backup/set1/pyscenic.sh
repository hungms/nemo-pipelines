#!/bin/bash -l

#SBATCH --job-name=pyscenic_set1
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export logpath=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/
mkdir -p ${logpath}pyscenic
rm -r /camp/home/hungm/template_pips/pyscenic/set1/slurm-*
cp -r /camp/home/hungm/template_pips/pyscenic/set1/ ${logpath}/pyscenic/
exec > ${logpath}pyscenic/set1/pyscenic.log 2>&1

conda activate /camp/home/hungm/.conda/envs/pyscenic
python ${logpath}pyscenic/set1/pyscenic.py
