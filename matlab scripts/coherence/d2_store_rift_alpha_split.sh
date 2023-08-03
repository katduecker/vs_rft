#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:00:0
#SBATCH --mem 100G
#SBATCH --qos bbdefault
#SBATCH --account=jenseno-visual-search-rft

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d2_store_rift_alpha_split([-1 0],0,1), quit"