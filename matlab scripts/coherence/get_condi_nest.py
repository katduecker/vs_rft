# -*- coding: utf-8 -*-
"""
extract condition values from CSV file to allow dynamic batch submission on bluebear
for condition specs with more than 1 entry
outputs: 
    - empty space separated string containing:
        - condi: condition specs
        - fw: width of bandpass filter in Hz
        - rej_sac: indicate if trials containing eye movement should be rejected
    
"""


import sys
import csv
from itertools import islice

def getcondnest(filepth,jobid):
    
    condi1 = []
    condi2 = []
    fw = []
    rej_sac = []
    
    with open(filepth) as f:
        reader = csv.DictReader(f, delimiter=',')
        for row in reader:
           # print(row)
            condi1.append(row["condi1"])
            condi2.append(row["condi2"])
            fw.append(row["fw"])
            rej_sac.append(row["rej_sac"])
            
       # fw = list(map(int,fw))
        #rej_sac = list(map(int,rej_sac))
        
        o = [condi1[jobid],condi2[jobid],fw[jobid],rej_sac[jobid]]
        return ' '.join(o)
        #return fw[jobid]
    
# this part is needed to parse (=split) input into separate arguments
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Get values from CSV file to pass to matlab')
    parser.add_argument('filepth', help="CSV file to read")
    parser.add_argument('jobid', type=int, help="index of job")
    args = parser.parse_args()

# print -> return to bash
    print(getcondnest(args.filepth, args.jobid))