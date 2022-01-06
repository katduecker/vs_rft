#!/usr/bin/env python
"""
Submit a separate Slurm job for each subject in the analysis

Submit some matlab code to the cluster for each subject.
In the code, replace the subject number with `<N>`
All other arguments get passed to sbatch_submit.py

Example
    # Double/single quotes must be preserved!
    python slurm_party.py "rs_preproc(<N>, 'trial')" -t 5:0 -m 20G

"""

import sys
import csv
import os

exp_dir = '/rds/projects/2018/jenseno-entrainment/subjects/Batch_3/'
fname = 'subject_info.csv'
number_tag = '<N>'

def command_template():
    ''' Return the template command with a variable for subject number
    '''
    cmd = 'ent_sbatch_submit.py -i "{}"'.format(sys.argv[1])
    n_args = len(sys.argv)
    for i_arg in range(2, n_args): # ignore script name (0) and matlab code (1)
        cmd = cmd + ' ' + sys.argv[i_arg]
    print cmd
    return cmd

def main():
    cmd = command_template()

    # Load the subject info
    with open(exp_dir + fname) as f:
        csv_reader = csv.DictReader(f)
        subject_info = [line for line in csv_reader]

    # Run the command for each subject
    for subj in subject_info:
        if subj['exclude'] == '1':
            continue
        else:
            c = cmd.replace(number_tag, subj['n'])
            print c
            os.system(c)

if __name__ == '__main__':
    main()
