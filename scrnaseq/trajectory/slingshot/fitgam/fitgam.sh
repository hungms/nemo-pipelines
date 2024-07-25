#!/bin/bash -l

#SBATCH --job-name=evaluatek
#SBATCH --partition=cpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000

export log_path=/camp/home/hungm/working/Matthew/project/anqi/BMPC_timestamped/logs/
cp -r /camp/home/hungm/template_pips/slingshot ${log_path}
exec > ${log_path}/slingshot/fitgam.log 2>&1

conda activate /camp/home/hungm/.conda/envs/slingshot
Rscript /camp/home/hungm/template_pips/slingshot/evaluatek.sh
