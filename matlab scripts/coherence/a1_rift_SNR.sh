#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 6:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 2
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a1_rift_SNR(${SLURM_ARRAY_TASK_ID},50:75,5,{'but','twopass'},'MEGGRAD'), quit"
