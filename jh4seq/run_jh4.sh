#!/bin/bash
#SBATCH --job-name=jh4
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:0
#SBATCH --mem=50G
#SBATCH --partition=ncpu
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=matthew.hung@crick.ac.uk

### edit this ###
export accession=oscar/OA_DN19023
export PRJ=/camp/home/hungm/working/Matthew/project/${accession}
#################

source /camp/home/hungm/pipelines/piplog.sh
ml Perl/5.26.1-GCCcore-6.4.0
export PATH=$PATH:/nemo/lab/caladod/working/Matthew/library/perl5/lib/perl5
export PERL5LIB=/nemo/lab/caladod/working/Matthew/library/perl5/lib/perl5:$PERL5LIB
perl -I/nemo/lab/caladod/working/Matthew/library/perl5 ${PRJ}/logs/jh4seq/run_jh4.pl
