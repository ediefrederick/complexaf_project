#!/bin/bash

# Arguments
FASTA_FILE=$1
ALPHAFOLD_LOG_FILE=$2


# Get the basename of the fasta file (e.g., file.fasta -> file)
fasta_basename=$(basename "$FASTA_FILE" .fasta)

# Define tmux session name based on the fasta file basename
session_name="alphafold_$fasta_basename"

# Start a new tmux session to run AlphaFold
tmux new -d -s "$session_name" "bash -c './run_alphafold.sh "$FASTA_FILE" "$ALPHAFOLD_LOG_FILE"'" 

