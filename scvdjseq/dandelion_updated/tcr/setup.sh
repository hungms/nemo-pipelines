export accession=/lingling/LZ_SC23355/
export cellranger_outdir=/camp/home/hungm/scratch/hungm/${accession}/cellranger/
export dandelion_dir=/camp/home/hungm/scratch/hungm/${accession}/dandelion/tcr/

mkdir -p ${dandelion_dir}

for sample in "$cellranger_outdir"/*;
do 
   if [ -d "${sample}" ]; then
      i=$(basename "${sample}")
      mkdir -p ${dandelion_dir}${i}

############# all contigs #############
      cp ${cellranger_outdir}${i}/outs/multi/vdj_t/all_contig.fasta ${dandelion_dir}${i}/
      cp ${cellranger_outdir}${i}/outs/multi/vdj_t/all_contig_annotations.csv ${dandelion_dir}${i}/
#######################################

############# filtered contigs #############
      #cp ${cellranger_outdir}${i}/outs/per_sample_outs/${i}/vdj_t/filtered_contig.fasta ${dandelion_dir}${i}/
      #cp ${cellranger_outdir}${i}/outs/per_sample_outs/${i}/vdj_t/filtered_contig_annotations.csv ${dandelion_dir}${i}/
############################################

   fi
done
