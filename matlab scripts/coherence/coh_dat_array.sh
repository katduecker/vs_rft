#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 3:0:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --array 2

set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {5..36}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d_extract_data_coh(${SLURM_ARRAY_TASK_ID},$i), quit"
done