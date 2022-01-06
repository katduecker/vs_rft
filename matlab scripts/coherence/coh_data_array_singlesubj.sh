#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 10:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 5-36

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d_extract_data_coh(3,${SLURM_ARRAY_TASK_ID}), quit"
