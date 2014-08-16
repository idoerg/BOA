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
import faa
import gff
import cPickle


if __name__=="__main__":
    
    
    parser = argparse.ArgumentParser(description=\
        'Extracts orfs and creates fasta files')
    parser.add_argument('--gff', type=str, required=False,
                        help='Input gff file')
    parser.add_argument('--fasta', type=str, required=False,
                        help='Input fasta file')
    parser.add_argument('--faidx', type=str, required=False,
                        help='Input fasta file index')
    parser.add_argument('--faa', type=str, required=False,
                        help='Input fasta file for proteins')
    parser.add_argument('--pickle', type=str, required=False,
                        help='Pickle file containing all clusters')
    parser.add_argument('--out', type=str, required=False,
                        help='Output file for translated orfs')
    args = parser.parse_args()
    queries = []
    if not os.path.exists("faa.pickle"):
        print "Creating pickle"
        faaindex = faa.FAA(args.faa)
        faaindex.index()
        print len(faaindex.indexer)
        assert 'AAD07798.1' in faaindex.indexer,"len=%d"%len(faaindex.indexer) 
        faaindex.dump("faa.pickle")
    else:
        print "Loading pickle"
        faaindex = faa.FAA(args.faa)
        faaindex.load("faa.pickle")
        
    gff = gff.GFF(args.gff,fasta_file=args.fasta,fasta_index=args.faidx)                
    clusters,predclusters = cPickle.load(open(args.pickle,'rb'))
    for cluster in clusters:
        for node in cluster:
            acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description = node.split('|')
            function = clrname.split('.')[0]
            if function=='toxin':
                queries.append((acc,clrname,full_evalue,hmm_st,hmm_end,env_st,env_end,description))        
        
        
    gff.translate_orfs(queries,faaindex,args.out)




        
        