#!/bin/bash -l

#SBATCH --job-name=nmf
#SBATCH --partition=ncpu
#SBATCH --ntasks=32
#SBATCH --time=3-00:00
#SBATCH --mem=250000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### change data path ###
export accession=oscar/OA_SC23423
export sample=20240711_integrated
export outdir=/camp/home/hungm/working/Matthew/project/${accession}/output/nmf/
export input_h5ad=/camp/home/hungm/working/Matthew/project/${accession}/mudata/merged/SC23423_merged_rchop_placebo_tumour_integrated_raw.h5ad
#export hvf=/camp/home/hungm/working/Matthew/project/anqi/AX_SC22272/seurat/SC22272_variablefeatures.txt
#export input_tpm=/camp/home/hungm/working/Matthew/project/anqi/AX_SC22272/anndata/integrated_data.txt
########################
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}/
source /camp/home/hungm/pipelines/piplog.sh
conda activate /camp/home/hungm/.conda/envs/nmf

cnmf prepare --output-dir ${outdir} --name ${sample} -k 2 3 4 5 6 7 8 9 10 11 12 13 --n-iter 100 --seed 14 --numgenes 3000  -c ${input_h5ad} #--tpm ${input_tpm} --genes-file ${hvf} #-c ${input_h5ad}
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 0 --total-workers 1
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 1 --total-workers 2
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 2 --total-workers 3
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 3 --total-workers 4
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 4 --total-workers 5
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 5 --total-workers 6
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 6 --total-workers 7
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 7 --total-workers 8
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 8 --total-workers 9
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 9 --total-workers 10
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 10 --total-workers 11
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 11 --total-workers 12
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 12 --total-workers 13
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 13 --total-workers 14
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 14 --total-workers 15
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 15 --total-workers 16
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 16 --total-workers 17
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 17 --total-workers 18
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 18 --total-workers 19
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 19 --total-workers 20
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 20 --total-workers 21
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 21 --total-workers 22
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 22 --total-workers 23
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 23 --total-workers 24
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 24 --total-workers 25
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 25 --total-workers 26
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 26 --total-workers 27
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 27 --total-workers 28
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 28 --total-workers 29
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 29 --total-workers 30
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 30 --total-workers 31
cnmf factorize --output-dir ${outdir} --name ${sample} --worker-index 31 --total-workers 32
cnmf combine --output-dir ${outdir} --name ${sample}
cnmf k_selection_plot --output-dir ${outdir} --name ${sample}

### run separate
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 3 --local-density-threshold 0.01 --show-clustering
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 4 --local-density-threshold 0.01 --show-clustering
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 5 --local-density-threshold 0.01 --show-clustering
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 6 --local-density-threshold 0.01 --show-clustering
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 7 --local-density-threshold 0.01 --show-clustering
#cnmf consensus --output-dir ${outdir} --name ${sample} --components 8 --local-density-threshold 0.01 --show-clustering
