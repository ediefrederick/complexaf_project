#!/bin/bash

# Set paths
FASTA_FOLDER="/mnt/users/edie01/AF/complex_fasta"
MAIN_LOG_FILE="/mnt/users/edie01/AF/alphafold_logs/main.log"


# Redirect stdout and stderr to the main log file (overwrite existing file)
exec >> "$MAIN_LOG_FILE" 2>&1
echo Starting!

# Loop through each fasta file in the folder
for fasta_file in "$FASTA_FOLDER"/*.fasta; do
  # Get the basename of the fasta file (e.g., file.fasta -> file)
  echo "Processing file: $fasta_file"
  fasta_basename=$(basename "$fasta_file" .fasta)

  # Define a separate log file for run_alphafold.sh
  ALPHAFOLD_LOG_FILE="/mnt/users/edie01/AF/alphafold_logs/${fasta_basename}_alphafold.log"

  # Run AlphaFold by calling run_alphafold.sh and pass the log file
  ./run_tmux.sh "$fasta_file" "$ALPHAFOLD_LOG_FILE"

  # Check if AlphaFold output directory is empty
#   if [ ! "$(ls -A "$OUTPUT_BASE_DIR")" ]; then
#     echo "AlphaFold run for $fasta_file failed. Output directory is empty."
#   fi
done
