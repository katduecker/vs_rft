#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 45:0
#SBATCH --qos bbdefault
#SBATCH --mem 80G

set -e

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_coh_hilb, quit"