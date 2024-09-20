#!/bin/env bash

#The first line (line 1) starting with HASH (#) is NECESSARY and tells the system that this is a bash shell script.
#It is possible (eg in the case a scheduler is available) to have more CONSECUTIVE HASHED (#) lines following (not used here).
#These lines will all be read when running the bash shell script.

#Once there is a line in the bash shell script not starting with a HASH (#) (line 2 in this sample script), from then on
#a HASH (#) is used to only start comment lines like this one. These lines will all be ignored when running the bash shell script.

#This sample script:
#1. Creates a new directory with an empty text file in it
#2. Runs a python script we already have in the same directory and print its result in standard output
#3. Redirects all standard outputs and errors to a custom log file. If the log file does not exits, it creates it.

#Redirect all "stderr" ("2") to "stdout" ("1") to "sample_script.log". create the log file if it does not exist, or append the redirected 
#output to the existing file. It will overwrite an existing file if changed to: exec > sample_script.log 2>&1
exec >> sample_script.log 2>&1

#Because of the exec redirection the output from ALL LINES until the end of the script (as well as any errors) will be written to sample_script.log
echo ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#We can use environmental variables as normal in a bash shell script. In the example bellow `date` is a system variable.
echo +++++++++++++++++++++++++	`date`	:	Starting Sample.sh script 	+++++++++++++++++++++++++

#Run more bash commands as you would use them during an interactive session at the Terminal
source /opt/anaconda3/bin/activate
conda  init
source venv/bin/activate

#Give variable names to standard input(s) passed when calling this script. We can pass any number of variables as standard input.
#Example: ./sample_script.sh some_name some_path
name=$1
path=$2

mkdir "$path"
touch "$path/$name.txt"
chmod 770 "$path/$name.txt"

#Call another script and write its output to a variable.
#It will return a "No such file or directory" error without a script2.py
new_variable=($(python3 script2.py))
echo $new_variable

echo +++++++++++++++++++++++++	`date`	:	Finished Sample.sh script 	+++++++++++++++++++++++++
