#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 3:00:0
#SBATCH --mem 150G
#SBATCH --qos bbdefault

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_fun_coh_cond(1,'6760',53:78,2), quit"