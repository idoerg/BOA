
import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess
import shutil
import time
from collections import defaultdict

class Cluster(object):
    def __init__(self):
        self.seqs = []
    def append(self,seq):
        self.seqs.append(seq)
    def __len__(self):
        return len(self.seqs)

class CDHit(object):
    def __init__(self,input_file,output_file,similarity):
        self.input = input_file
        self.output = output_file
        self.similarity = similarity
        self.cluster_out = "%s.clstr"%(output_file)
        self.cluster_counts = "%s.count"%(output_file)
        self.clusters = list()
    def __len__(self):
        return len(self.clusters)    
    #Cleans up all of the temporary files generated
    def cleanup(self):
        pass

    """
    Runs CD hit script
    """
    def run(self):
        cmd = "cdhit -i %s -o %s -d 150 -c %f" %(self.input,self.output,self.similarity)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()

    cluster_reg = re.compile(r"")


    """
    Parses clustering information from CD hit
    """
    def parseClusters(self):
        self.clusters = list()
        with open(self.cluster_out,'r') as handle:
            for ln in handle:
                ln = ln.rstrip()
                if ln[0]==">":
                    self.clusters.append( Cluster() )
                else:
                    self.clusters[-1].append(ln)
    def countOut(self):
        handle = open(self.cluster_counts,'w')
        for clr in self.clusters:
            handle.write("%d\n"%(len(clr)))
    def filterSize(self,n):
        self.clusters = [clr for clr in self.clusters if len(clr)>n]
        
def go(input_file,output_file,threshold):
    cdhitProc = CDHit(input_file,output_file,threshold)
    cdhitProc.run()
    cdhitProc.parseClusters()
    cdhitProc.countOut()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Generates counts for cluster')
    parser.add_argument(\
        '--input', type=str, required=False,
        help='Input fasta file')
    parser.add_argument(\
        '--output', type=str, required=False,
        help='Output cluster counts')
    parser.add_argument(\
        '--threshold', type=float, required=False, default=0.9,
        help='Similarity threshold score for clustering')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run unittests')
    args = parser.parse_args()

    if not args.test:
        go(args.input,args.output,args.threshold)
    else:
        del sys.argv[1:]
        import unittest
        import test
        class TestCluster(unittest.TestCase):
            def testcluster1(self):
                pass
        unittest.main()

