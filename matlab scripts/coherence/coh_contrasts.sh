#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 6:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 3

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d2_coh_contrasts(${SLURM_ARRAY_TASK_ID},'6067'), quit"