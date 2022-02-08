#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 2:00:0
#SBATCH --mem 120G
#SBATCH --qos bbdefault
#SBATCH --array 2

set -eu

module purge; module load bluebear
module load MATLAB/2019b

for i in {0..3}; do
	condi=($(python3 get_condi_nest.py "coh_condi_nest.csv" $i))
	echo ${condi[0]} 
	echo ${condi[1]} 
	echo ${condi[2]} 
	echo ${condi[3]} 
	#matlab -nodisplay -r "run /rds/homes/d/dueckerk/startup.m, a_fun_coh_cond_new(${SLURM_ARRAY_TASK_ID},{'${condi[0]}','${condi[1]}'},53:78,${condi[2]},${condi[3]}), quit"
done