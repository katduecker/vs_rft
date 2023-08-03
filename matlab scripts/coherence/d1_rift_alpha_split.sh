#!/bin/bash
#SBATCH --ntasks 4                               
#SBATCH --time 2:00:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-33
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {1..4}; do
	matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d1_rift_alpha_split(${SLURM_ARRAY_TASK_ID},$i,5,[-1 0],0), quit"
done