#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:0
#SBATCH --mem=8G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

####### edit here #########
export accession=ellie/EA_RN24134
###########################

mkdir -p /camp/home/hungm/scratch/hungm/${accession}/
export PRJ=/camp/home/hungm/scratch/hungm/${accession}/
source /camp/home/hungm/pipelines/piplog.sh

# many of the changes you need to make to this file are the same as those 
# changes made to the '1.fetchngs-job-template-rnaseq.sh' file. Please refer to
# the comments in '1.fetchngs-job-template-rnaseq.sh' for tips on modifying
# this file.

# Load software modules - check for newer versions using 'ml spider Nextflow'
ml purge
ml Nextflow/23.10.0
ml Singularity/3.6.4
ml CAMP_proxy
export _JAVA_OPTIONS=-Djava.io.tmpdir=/flask/scratch/caladod/hungm/tmp

# Create a 'nextflow' folder in your working area. Run 'cd nextflow' and 'pwd -P'
# Copy the result to replace the directory path below
export NXF_HOME=/flask/scratch/caladod/hungm/nextflow


# NXF_WORK will be created in the directory you execute this script from
export NXF_WORK=/flask/scratch/caladod/hungm/work/


# You can leave this line as is - we are using the BABS cache of singularity images
export NXF_SINGULARITY_LIBRARYDIR=/flask/apps/containers/all-singularity-images/

# You should create a folder in your working area to store any additional images
# that are not found in the BABS library folder above
export NXF_SINGULARITY_CACHEDIR=/flask/scratch/caladod/hungm/singularity_image/

# run the pipeline
nextflow run nf-core/rnaseq \
    --input /flask/scratch/caladod/hungm/${accession}/input/paired_end.csv  \
    --outdir /flask/scratch/caladod/hungm/${accession}/output/ \
    --aligner star_salmon \
    --fasta /flask/reference/Genomics/babs_nfcore/genomes/mus_musculus/ensembl/GRCm39/genome/genome.fa \
    --gtf /flask/reference/Genomics/babs_nfcore/genomes/mus_musculus/ensembl/GRCm39/annotation/release-105/gtf/annotation.gtf \
    --email matthew.hung@crick.ac.uk \
    --skip_biotype_qc \
    --save_reference \
    -profile crick \
    -r 3.13.2 \
    -resume

## change genome reference
#    --fasta /flask/reference/Genomics/babs_nfcore/genomes/mus_musculus/ensembl/GRCm39/genome/genome.fa \
#    --gtf /flask/reference/Genomics/babs_nfcore/genomes/mus_musculus/ensembl/GRCm39/annotation/release-105/gtf/annotation.gtf \
#    --fasta /flask/reference/Genomics/babs_nfcore/genomes/homo_sapiens/ensembl/GRCh38/genome/genome.fa \
#    --gtf /flask/reference/Genomics/babs_nfcore/genomes/homo_sapiens/ensembl/GRCh38/annotation/release-105/gtf/annotation.gtf \


# Once to have made your changes to this file, save it and submit the job with:
#   sbatch 2.rnaseq-job-template.sh 

