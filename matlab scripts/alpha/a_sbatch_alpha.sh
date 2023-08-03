#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 1:00:0
#SBATCH --mem 60G
#SBATCH --qos bbdefault
#SBATCH --array 1-33
#SBATCH --account=jenseno-entrainment
set -eu

module purge; module load bluebear
module load MATLAB/2019b

matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_fun_alpha_pow(${SLURM_ARRAY_TASK_ID},1,4:30,[-2.25 1.75],0,0.5), quit"

