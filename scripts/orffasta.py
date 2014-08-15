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
import cPickle


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Extracts orfs and creates fasta files')
    parser.add_argument('--fasta', type=str, required=False,
                        help='Input fasta file')
    parser.add_argument('--faidx', type=str, required=False,
                        help='Input fasta index')
    parser.add_argument('--pickle', type=str, required=False,
                        help='Pickle file containing all clusters')
    
    args = parser.parse_args()
    indexer = fasta.Indexer(args.fasta,args.faidx)
    clusters,predclusters = cPickle.load(open(args.pickle,'rb'))
    indexer.load()
    for cluster in clusters:
        for node in cluster:
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description = node.split('|')
            function = clrname.split('.')[0]
            if function=="toxin":
                acc = "_".join(acc.split('_')[:-1])
                st,end = map(int,[env_st,env_end])
                seq = indexer.fetch(acc,st,end)
                print ">acc=%s|start=%d|end=%d|%s\n%s"%(acc,st,end,description,fasta.format(seq))
    
        
        