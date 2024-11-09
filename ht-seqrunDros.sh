BAM_DIR="/home/gerardo/DrosPepper"
# Define the GTF file
GTF_FILE="/home/gerardo/Drosophila_melanogaster.BDGP6.46.113.gtf/Drosophila_melanogaster.BDGP6.46.113.gtf"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
  # Define the output file name by replacing .bam with _counts.txt
  OUTPUT_FILE="${BAM_FILE%.bam}_counts.txt"
  
  # Print the BAM file being processed
  echo "Processing $BAM_FILE..."
  
  # Run htseq-count
  htseq-count -f bam -r pos -s no -t gene -i gene_id "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_FILE"
  
  # Check if htseq-count was successful
  if [ $? -eq 0 ]; then
    echo "Counts for $BAM_FILE written to $OUTPUT_FILE"
  else
    echo "Error processing $BAM_FILE"
  fi
done

echo "htseq-count analysis completed for all samples."