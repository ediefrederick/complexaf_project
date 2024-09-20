#!/bin/bash

# Define variables
FASTA_FILE="/mnt/users/edie01/AF/complex_fasta/T1106.fasta"
OUTPUT_DIR="/mnt/users/edie01/AF/complex_output"
LOG_FILE="/mnt/users/edie01/AF/alphafold_logs/alphafold_simple.log"
SESSION_NAME="alphafold_session"

# Redirect stdout and stderr to the specified log file (overwrite existing file)
exec >> "$LOG_FILE" 2>&1

# Start a new tmux session to run AlphaFold
tmux new -d -s "$SESSION_NAME" "bash -c 'python3 /mnt/scratch/alphafold/docker/run_docker.py \
    --fasta_paths=\"$FASTA_FILE\" \
    --max_template_date=2022-01-01 \
    --data_dir=/mnt/scratch/alphafold_ref \
    --output_dir=\"$OUTPUT_DIR\"'"

# Output the tmux session status
echo "Tmux session '$SESSION_NAME' started with AlphaFold command."
