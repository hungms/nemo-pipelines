#!/bin/bash -l

#SBATCH --job-name=pyscenic
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export log_path=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/
cp -r /camp/home/hungm/template_pips/pyscenic ${log_path}
exec > ${log_path}/pyscenic/pyscenic.log 2>&1

conda activate /camp/home/hungm/.conda/envs/pyscenic
python ${log_path}/pyscenic/pyscenic.py
