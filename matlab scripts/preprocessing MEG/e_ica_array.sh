#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 30:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --array 1-48

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, e1_fun_ICA(${SLURM_ARRAY_TASK_ID},0), quit"
