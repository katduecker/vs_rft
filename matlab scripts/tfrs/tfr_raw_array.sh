#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:30:0
#SBATCH --mem 70G
#SBATCH --qos bbdefault
#SBATCH --array 3-8

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, c1_tfr_raw(${SLURM_ARRAY_TASK_ID}), quit"