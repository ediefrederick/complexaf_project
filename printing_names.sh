#!/bin/bash

# Set paths
FASTA_FOLDER="/mnt/users/edie01/AF/complex_fasta"
LOG_FILE="/mnt/users/edie01/alphafold_logs/sample_script.log"

# Redirect stdout and stderr to the log file (overwrite existing file)
exec > "$LOG_FILE" 2>&1

# Loop through each fasta file in the folder
for fasta_file in "$FASTA_FOLDER"/*.fasta; do
  # Print the name of each fasta file
  echo "Processing file: $fasta_file"
done

echo "File listing completed."
