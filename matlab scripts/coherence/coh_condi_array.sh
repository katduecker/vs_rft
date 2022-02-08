#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 3:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 0-4

set -eu

module purge; module load bluebear
module load MATLAB/2019b

condi=($(python3 get_condi.py "coh_condi_tfreq.csv" ${SLURM_ARRAY_TASK_ID}))

#echo ${condi[0]}
matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_fun_coh_cond_new(1,{'${condi[0]}'},53:78,${condi[1]},${condi[2]}), quit"