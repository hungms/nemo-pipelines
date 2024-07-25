#!/bin/bash -l

#SBATCH --job-name=pyscenic_set4
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export logpath=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/
mkdir -p ${logpath}pyscenic
rm -r /camp/home/hungm/template_pips/pyscenic/set4/slurm-*
cp -r /camp/home/hungm/template_pips/pyscenic/set4/ ${logpath}/pyscenic/
exec > ${logpath}pyscenic/set4/pyscenic.log 2>&1

conda activate /camp/home/hungm/.conda/envs/pyscenic
python ${logpath}pyscenic/set4/pyscenic.py
