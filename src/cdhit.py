
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

def addArgs(parser):
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file basename containing the cdhit information')
    parser.add_argument(\
        '--keep-tmp', action='store_const', const=True, default=False,
        help='Keeps temporary files such as blast database and blast output xml')
    parser.add_argument(\
        '--threshold', type=float, required=False, default=0.9,
        help='Similarity threshold score for clustering')

class CDHit(object):
    def __init__(self,input_file,output_file,similarity):
        self.input = input_file
        self.output = output_file
        self.similarity = similarity
        self.cluster_output = "%s.clstr"%(output_file)
    #Cleans up all of the temporary files generated
    def cleanup(self):
        pass
    """
    Runs CD hit script
    """
    def run(self):
        cmd = "cdhit -i %s -o %s -c %f" %(self.input,self.output,self.similarity)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()
        
    cluster_reg = re.compile(r"")
    def processCluster(self):
        pass
    """
    Parses clustering information from CD hit
    """
    def parseClusters(self):
        buf = queue()
        cluster = defaultdefault( list )
        
        with open(self.cluster_out,'r') as handle:
            for ln in handle:
                ln = ln.rstrip()
                if ln[0]!=">":
                    clusterNum = ln[1:]
                    cluster[clusterNum].append(ln)
        
        
        
        pass

if __name__=="__main__":
    import unittest
    import test
    class TestCluster(unittest.TestCase):
        def testcluster1(self):
            pass
    unittest.main()
