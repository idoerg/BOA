
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
import cPickle
import gzip
import copy
from collections import defaultdict

class Cluster(object):
    def __init__(self):
        self.seqs = []
        self.iterator = iter(self.seqs)
    def append(self,seq):
        self.seqs.append(seq)
    def __len__(self):
        return len(self.seqs)
    def __str__(self):
        return '\n'.join(self.seqs)
    def __iter__(self):
        return self
    
    def next(self):
        try:
            record = next(self.iterator)
            return record
        except StopIteration as s:
            raise StopIteration()
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
    def __str__(self):
        return '\n'.join(map(str,self.clusters))
    """Cleans up all of the temporary files generated"""
    def cleanup(self):
        pass
    """ Dump object into pickle file """
    def dump(self,outfile):
        cPickle.dump(self.clusters,gzip.GzipFile(outfile,'wb'))
        #cPickle.dump(self,open(outfile,'wb'))
    """ Load object from pickle file """
    def load(self,infile):
        self.clusters = cPickle.load(gzip.GzipFile(infile,'rb'))
        #self = cPickle.load(open(infile,'rb'))
        
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
        self.clusters = [ clr for clr in self.clusters if len(clr)>n ]
        
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
            def setUp(self):
                seqs = [">test1",
                        "ACGTAACGTAACGTACGTACGTACGTACGTACGT",
                        ">test2",
                        "ACGTAACGTAACGTACGGACGTACGTACGTACGT",
                        ">test3",
                        "TGACTGCTGACTGCTGACTGCTGACTGCTGGTGGGG",
                        ">test4",
                        "TGACTGCTGACTGCTGACTGCTGACTGCTGGGGGGG"
                        ]
                self.fasta = "test.fa"
                self.clusterfile = "clustered"
                self.clusterreps = "clustered.clstr"
                open(self.fasta,'w').write('\n'.join(seqs))
                self.zip = "test_serial.zip"
            def tearDown(self):
                os.remove(self.clusterfile)
                os.remove(self.clusterreps)
                os.remove(self.fasta)
                
            def testcluster1(self):
                cdhitproc = CDHit(self.fasta,self.clusterfile,0.8)
                cdhitproc.run()
                cdhitproc.parseClusters()
                print cdhitproc
                self.assertEquals(len(cdhitproc.clusters),2)
                
            def testSerialization(self):
                cdhitproc = CDHit(self.fasta,self.clusterfile,0.8)
                cdhitproc.run()
                original = CDHit(self.fasta,self.clusterfile,0.8)
                cdhitproc.parseClusters()
                cdhitproc.dump(self.zip)
                original.load(self.zip)
                self.assertEquals(len(original.clusters),len(cdhitproc.clusters))
                self.assertEquals(str(original),str(cdhitproc))
                os.remove(self.zip)
        unittest.main()

