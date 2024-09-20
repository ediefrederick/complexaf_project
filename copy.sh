#!/bin/bash

# Define variables
REMOTE_USER="edie01"            
REMOTE_HOST="129.85.62.21"          
REMOTE_FILE_PATH="/mnt/users/edie01/AF/complex_output/T1160"
LOCAL_DESTINATION="/Users/edie/Desktop/thesis"

# Use scp to copy the file
scp -r -P 2013  "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_FILE_PATH}" "${LOCAL_DESTINATION}"

