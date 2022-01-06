#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 30:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault

set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, d_extract_data_coh(2,5), quit"