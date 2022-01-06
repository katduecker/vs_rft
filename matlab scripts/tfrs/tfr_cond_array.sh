#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 3:0:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --array 1-8

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, c2_tfr_Tfreq(${SLURM_ARRAY_TASK_ID}), quit"