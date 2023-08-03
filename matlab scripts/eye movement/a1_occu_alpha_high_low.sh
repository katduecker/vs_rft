#!/bin/bash
#SBATCH --time 1:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft
#SBATCH --array 1-31

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a1_occu_alpha_high_low(${SLURM_ARRAY_TASK_ID},[0.25 0.5]), quit"