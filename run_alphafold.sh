#!/bin/bash

# Set paths
OUTPUT_DIR="/mnt/users/edie01/AF/complex_output_second"

# Arguments
FASTA_FILE=$1
ALPHAFOLD_LOG_FILE=$2

# Redirect stdout and stderr to the specified log file (overwrite existing file)
exec >> "$ALPHAFOLD_LOG_FILE" 2>&1

# Run AF
python3 /mnt/scratch/alphafold/docker/run_docker.py \
    --fasta_paths="$FASTA_FILE" \
    --max_template_date=2022-01-01 \
    --model_preset=multimer \
    --data_dir=/mnt/scratch/alphafold_ref \
    --output_dir="$OUTPUT_DIR"