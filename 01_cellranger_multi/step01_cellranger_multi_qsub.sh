#!/bin/bash

#----------- define global variables
base_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025

# path to the cell ranger script
cellranger_script=$base_dir/scripts/GB-SB-1467/01_cellranger_multi/cellranger_multi.sh

# output directory
outdir=$base_dir/results/01_cellranger_multi/

# create the results output directory
mkdir -p $outdir


#----------- run cell ranger multi script for all config files
# Find all CSV config files
for file in $base_dir/scripts/GB-SB-1467/01_cellranger_multi/*.csv; do
    # Extract the file name from the path
    file_name=$(basename "$file")
    
    # Remove the _config.csv extension to get the base name
    base_name="${file_name%.csv}"
    sample_out_id="${base_name%_config}"

    # print variables for logs
    echo "Processing config file: $file, Output name: $sample_out_id, Output path: $outdir"

    # run cellranger multi for current sample
    qsub $cellranger_script $sample_out_id $file $outdir
done


################### END ###################
