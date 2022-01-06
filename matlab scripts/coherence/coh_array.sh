#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 70:0
#SBATCH --mem 80G
#SBATCH --qos bbdefault
#SBATCH --array 2

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_coh_hilb(${SLURM_ARRAY_TASK_ID}), quit"