#!/usr/bin/env bash

data_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
script_dir=$data_dir/scripts/GB-SB-1467/05_cluster_resolution_optimization

# Define the list of PCs
pcs_list=(30 20 15)

# Iterate over each value in pcs_list
for npcs in "${pcs_list[@]}"
do
    # Submit the job using qsub
    qsub $script_dir/qsub_05b_cluster_resolution_optimization_per_sample.sh "$npcs"
done

echo "Script completed!"

######################## END ########################