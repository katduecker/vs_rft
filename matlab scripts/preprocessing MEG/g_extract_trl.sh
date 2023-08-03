#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 20:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-34
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, g_split_trl(${SLURM_ARRAY_TASK_ID},0.05), quit"
