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
import gff

class FAA(object):
    def __init__(self,faafile,gfffile):
        self.faafile = faafile
        self.gfffile = gfffile
        gff_index = 0
    def index(self):
        for record in SeqIO.parse(open(self.faafile,'r'),'fasta'):
            #print record
            #print record.id
            #print record.id.split('|')
            gi,num,dbtype,protid,_ = record.id.split('|')
            #print gi,num,dbtype,protid
            self.indexer[protid] = record
        
    def dump(self,outfile):
        cPickle.dump((self.indexer,self.faafile),open(outfile,'wb'))
    def load(self,outfile):
        self.indexer,self.faafile = cPickle.load(open(outfile,'rb'))
    def __getitem__(self,protid):
        return self.indexer[protid]
    

    
def reformat(orgfaa,faafile):
    outhandle = open(faafile,'w')
    for record in SeqIO.parse(open(orgfaa,'r'),'fasta'):
        gi,num,dbtype,protid,_ = record.id.split("|")
        outhandle.write(">%s\n%s\n"%(protid,fasta.format(str(record.seq))))
    outhandle.close()
    
    
    

    
    
    