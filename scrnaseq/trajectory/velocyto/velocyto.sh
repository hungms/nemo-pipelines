#!/bin/bash -l

#SBATCH --job-name=velocyto
#SBATCH --partition=cpu
#SBATCH --ntasks=3
#SBATCH --time=3-00:00
#SBATCH --mem=100000

export accession=GSE189775
export sample=GSE189775
export inputdir=/camp/home/hungm/scratch/hungm/${accession}/cellranger/${sample}/outs/per_sample_outs/${sample}/count/
export outputdir=/camp/home/hungm/scratch/hungm/${accession}/velocyto/${sample}

export bam=${inputdir}sample_alignments.bam
export gtf=/camp/home/hungm/scratch/hungm/reference/alignment/velocyto/GRCh38_rmsk.gtf
export rna=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf

rm -r /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/slurm*
exec > /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${sample}/velocyto.log 2>&1

conda activate /camp/home/hungm/.conda/envs/velocyto
cd ${outdir}
velocyto run -m ${gtf} ${bam} ${rna}
mv ${inputdir}velocyto ${outputdir}


## change mouse to human
#export gtf=/camp/home/hungm/scratch/hungm/reference/alignment/velocyto/mm10_rmsk.gtf
#export rna=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-mm10-2020-A/genes/genes.gtf

#export gtf=/camp/home/hungm/scratch/hungm/reference/alignment/velocyto/GRCh38_rmsk.gtf
#export rna=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf
