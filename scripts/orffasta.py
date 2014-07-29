"""
Creates fasta file from orfs
"""

import argparse
import os
import sys
import site
import re
import numpy as np
import numpy.random

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base_path="%s/src"%base_path
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
    
import fasta

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Extracts orfs and creates fasta files')
    parser.add_argument(\
                        '--fasta', type=str, required=False,
                        help='Input fasta file')
    parser.add_argument(\
                        '--faidx', type=str, required=False,
                        help='Input fasta index')
    args = parser.parse_args()
    indexer = fasta.Indexer(args.fasta,args.faidx)
    indexer.load()
    for ln in sys.stdin:
        ln = ln.rstrip()
        toks = ln.split('|')
        six_acc = toks[0]
        acc = "_".join(six_acc.split('_')[:-1])
        st,end = map(int,toks[5:7])
        seq = indexer.fetch(acc,st,end)
        print ">%s\n%s"%(ln,fasta.format(seq))
        
        