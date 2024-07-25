export accession=/tony/TC_SC24124

mkdir -p /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger

while IFS= read -r id;
do mkdir -p /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}
   cp -r /camp/home/hungm/pipelines/scrnaseq/cellranger/multi.sh /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/
   cp -r /camp/home/hungm/pipelines/scrnaseq/cellranger/config.csv /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/
   cp -r /camp/home/hungm/pipelines/scrnaseq/cellranger/batch* /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/

   sed -i "11s|.*|export accession=$accession|" /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/multi.sh
   sed -i "12s|.*|export sample=$id|" /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/multi.sh

   #csv1="/camp/home/hungm/scratch/hungm/$accession/logs/fastq_path/${id}_fastq.csv"
   #csv1="/camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/fastq_path/${id}_path.csv"
   #sed -i '37r '"$csv1"'' /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/config.csv

   #replacement="reference,/camp/home/hungm/scratch/hungm/$accession/input/${id}_feature_reference.csv"
   #sed -i "26s|.*|$replacement|" /camp/home/hungm/scratch/hungm/${accession}/logs/cellranger/${id}/config.csv;
done < /camp/home/hungm/pipelines/scrnaseq/cellranger/batch_id.txt

#########################
### batch submit jobs ###
#########################
#while IFS= read -r id;
#do sbatch ${id}/multi.sh;
#done < /camp/home/hungm/pipelines/scrnaseq/cellranger/batch_id.txt

#while IFS= read -r id;
#do rm -r /camp/home/hungm/scratch/hungm/${accession}/cellranger/${id};
#done < /camp/home/hungm/pipelines/scrnaseq/cellranger/batch_id.txt
