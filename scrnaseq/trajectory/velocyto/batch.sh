export accession=GSE230705

mkdir -p /camp/home/hungm/scratch/hungm/${accession}/logs/
mkdir -p /camp/home/hungm/scratch/hungm/${accession}/velocyto
mkdir -p /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto

while IFS= read -r id;
do ### make log file
   mkdir -p /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${id}
   cp -r /camp/home/hungm/template_pips/trajectory/velocyto/* /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${id}/
   rm -r /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${id}/batch*
   cp -r /camp/home/hungm/template_pips/trajectory/velocyto/batch* /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/

   ### edit script
   sed -i "9s/.*/export accession=$accession/" /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${id}/velocyto.sh
   sed -i "10s/.*/export sample=$id/" /camp/home/hungm/scratch/hungm/${accession}/logs/velocyto/${id}/velocyto.sh;
done < /camp/home/hungm/template_pips/trajectory/velocyto/batch_id.txt

#########################
### batch submit jobs ###
#########################
#while IFS= read -r id;
#do sbatch ${id}/velocyto.sh;
#done < /camp/home/hungm/template_pips/trajectory/velocyto/batch_id.txt
