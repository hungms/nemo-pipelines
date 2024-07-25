for file in *_1.fastq.gz; do
  newname=$(echo "$file" | sed 's/_1.fastq.gz/_S1_R1_001.fastq.gz/')
  mv "$file" "$newname"
done

for file in *_2.fastq.gz; do
  newname=$(echo "$file" | sed 's/_2.fastq.gz/_S1_R2_001.fastq.gz/')
  mv "$file" "$newname"
done

for file in *_3.fastq.gz; do
  newname=$(echo "$file" | sed 's/_3.fastq.gz/_S1_I1_001.fastq.gz/')
  mv "$file" "$newname"
done
