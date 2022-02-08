#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 2:30:0
#SBATCH --mem 80G
#SBATCH --qos bbdefault
#SBATCH --array 1-32
#SBATCH --account=jenseno-entrainment


set -eu

module purge; module load bluebear
module load MATLAB/2019b


condi=($(python3 get_condi_nest.py "coh_condi_nest.csv" 3))
#echo ${condi[0]} 
matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_fun_coh_cond_new(${SLURM_ARRAY_TASK_ID},{'${condi[0]}','${condi[1]}'},53:78,${condi[2]},${condi[3]}), quit"
