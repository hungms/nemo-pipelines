#!/bin/bash -l

#SBATCH --job-name=pyscenic_set3
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export logpath=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/
mkdir -p ${logpath}pyscenic
rm -r /camp/home/hungm/template_pips/pyscenic/set3/slurm-*
cp -r /camp/home/hungm/template_pips/pyscenic/set3/ ${logpath}/pyscenic/
exec > ${logpath}pyscenic/set3/pyscenic.log 2>&1

conda activate /camp/home/hungm/.conda/envs/pyscenic
python ${logpath}pyscenic/set3/pyscenic.py
