"""
Data structure that stores faa entries
"""

import os,sys,site
import numpy
base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
import fasta
import argparse
import hmmer
from Bio.Seq import Seq
from Bio import SeqIO
import re
from collections import defaultdict
from bx.intervals import *
import cPickle
class FAA(object):
    def __init__(self,faafile):
        self.faafile = faafile
        self.indexer = dict()
    def index(self):
        for record in SeqIO.parse(open(self.faafile,'r'),'fasta'):
            print record
            print record.id
            print record.id.split('|')
            gi,num,protid,dbtype,_ = record.id.split('|')
            self.indexer[protid] = record
    def dump(self,outfile):
        cPickle.dump((self.indexer,self.faafile),open(outfile,'w'))
    def load(self):
        self.indexer,self.faafile = cPickle.load(open(outfile,'r'))
    def __getitem__(self,protid):
        return self.indexer[protid]
    
