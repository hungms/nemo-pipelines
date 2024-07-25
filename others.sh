cat integrated_data.csv | sed 's/,/\t/g' | sed 's/\"//g' | awk '
BEGIN { OFS="\t" }
{
    for (i=1; i<=NF; i++) {
        a[i,NR] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) { 
        if (j == 1) {
            printf "%s", "";  # Print an empty string for the first row
        } else {
            printf "%s%s", a[j,1], OFS;  # Print the original column names
        }
        for(i=2; i<=NR; i++) { 
            printf "%s%s", a[j,i], (i<NR ? OFS : "\n");
        }
    }
}
' > integrated_data.txt

cat rna_counts.csv | sed 's/,/\t/g' | sed 's/\"//g' | awk '
BEGIN { OFS="\t" }
{
    for (i=1; i<=NF; i++) {
        a[i,NR] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) { 
        if (j == 1) {
            printf "%s", "";  # Print an empty string for the first row
        } else {
            printf "%s%s", a[j,1], OFS;  # Print the original column names
        }
        for(i=2; i<=NR; i++) { 
            printf "%s%s", a[j,i], (i<NR ? OFS : "\n");
        }
    }
}
' > rna_counts.txt


cat integrated_data.csv | sed 's/,/\t/g' | sed 's/\"//g' > integrated_data.txt
cat rna_counts.csv | sed 's/,/\t/g' | sed 's/\"//g' > rna_counts.txt















