submit jobs using
python_load='source /rds/homes/d/dueckerk/bin/load_python.sh'

sbatch_submit.py -i "python test.py" -s "$python_load" -t 2:0 -m 5G -a "jenseno-visual-search-rft"

sbatch_submit.py -i "python run_bio_mdl.py <delta> <p> <k>" -s "$python_load" -t 4:30:0 -m 80G -a "jenseno-visual-search-rft"

