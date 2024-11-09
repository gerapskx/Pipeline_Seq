#!/bin/bash

# Define the genome directory
genomedir="/home/gerardo/DrosDirectory_october"

# Get a list of unique sample prefixes

mylist="$(ls *.fastq.gz | perl -pe 's/^(.*?_SRR\d{8})_pass_\d_pair\.fastq\.gz$/$1/' | uniq)"

# Loop through each sample prefix
for prefix in $mylist; do
  pairedend1=`ls ${prefix}_pass_1_pair.fastq.gz`
  pairedend2=`ls ${prefix}_pass_2_pair.fastq.gz`

if [ -f "$pairedend1" ] && [ -f "$pairedend2" ]; then
    echo "Aligning sample: $prefix"
else
    echo "Warning: One or both files for sample $prefix do not exist."
    echo "Missing files: $pairedend1, $pairedend2"
    continue  # Skip to the next sample
fi

    # Align reads using STAR
    STAR \
        --runThreadN 8 \
        --genomeDir "$genomedir"  \
        --readFilesIn $pairedend1 $pairedend2 \
        --readFilesCommand zcat \
        --outFileNamePrefix "$prefix" \
        --outSAMtype BAM SortedByCoordinate \

echo 
'Alignment completed for:' $prefix 

done
