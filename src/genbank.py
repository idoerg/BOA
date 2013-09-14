"""
Wrapper class for genbank files
"""

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

import sys
import os
import argparse

class GenBank(object):
    def __init__(self,fname):
        self.fname = fname
        self.genes = dict()

        seq_record = SeqIO.parse(open(fname), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try:
                    gene = feature.qualifiers["gene"][0]
                    self.genes[ gene ]=feature
                except KeyError:
                    continue
    def findGene(self,gname):
        try:
            return self.genes[gname]
        except KeyError:
            print >> sys.stderr,"No such gene"

if __name__=="__main__":
    import unittest



    unittest.main()
