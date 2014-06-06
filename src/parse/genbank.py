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
    def __init__(self,fname,dicttype="gene"):
        self.fname = fname
        if dicttype=="protein":
            self.dict = self.buildProteinDictionary(fname)
        else:
            self.dict = self.buildGeneDictionary(fname)


    def buildProteinDictionary(self,fname):
        proteins = dict()
        print fname
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try:
                    protein = feature.qualifiers["protein_id"][0]
                    proteins[ protein ]=feature
                except KeyError:
                    continue
        return proteins

    def buildGeneDictionary(self,fname):
        genes = dict()
        seq_record = SeqIO.parse(open(fname), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try:
                    gene = feature.qualifiers["gene"][0]
                    genes[ gene ]=feature
                except KeyError:
                    continue
        return genes

    def findGene(self,gname):
        try:
            return self.dict[gname]
        except KeyError:
            print >> sys.stderr,"No such gene"

    def findProtein(self,pname):
        try:
            return self.dict[pname]
        except KeyError:
            print >> sys.stderr,"No such protein"

if __name__=="__main__":
    import unittest



    unittest.main()
