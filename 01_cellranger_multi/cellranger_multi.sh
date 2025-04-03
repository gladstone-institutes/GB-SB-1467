#!/bin/bash
#$ -cwd
#$ -o /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -e /gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025/tmp/
#$ -pe smp 4
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=20:00:00
#$ -j yes

###############################################################################
## Project ID: GB-SB-1467
## Authors: Ayushi Agrawal
##
## Script Goal: run cellranger multi
##
## Expected command line arguments:
##    args[1]: Name of the directory containing all cellranger-generated outputs
##    args[2]: Path of the config file
##    args[3]: Path of the output directory
##
## Example running this script on wynton:
##    qsub cellranger_multi.sh \
##	  gem_well_poll1 \
##	  config.csv \
##	  results/01_cellranger_multi
###############################################################################

#setup paths
base_dir=/gladstone/jain/boinformatics-collaboration/sb-1467-skyler-blume-isha-jain-snrnaseq-mm10-mar-2025
container_dir=$base_dir/assets
export APPTAINER_BINDPATH="$base_dir"

#assign the argument values to variables
sample_out_id=$1
csv_path=$2
outdir=$3

#chnage working directory to output directory
cd $outdir

#run cellranger count
apptainer exec $container_dir/cellranger_9.0.1.sif cellranger multi \
	--id=$sample_out_id \
	--csv=$csv_path


# End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"


## END 
