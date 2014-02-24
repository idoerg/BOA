"""
Finds intergenic regions in a bacterial genome
"""

import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone

from collections import defaultdict

import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess
import intervals

import genbank
import blast
import intergene

loc_reg = re.compile("(\S\S_\d+)-ign-\d+:(\d+)-(\d+)(\S)")

class IntergeneHandler:
    def __init__(self,intergene_file):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.intergene_file = intergene_file
        
        
    """
    Retrieves intervals that can be used to determine if bacteriocin overlap any intergenic regions
    """
    def getIntervals(self):
        self.intergeneDict = defaultdict(list)
        for record in SeqIO.parse(open(self.intergene_file,'r'),"fasta"):
            match = loc_reg.findall(record.id)[0]
            accession,start,end,strand = match
            self.intergeneDict[(accession,strand)].append( (int(start),int(end),strand) )            
    """
    Returns true if bacteriocin interval overlaps intergene region
    Otherwise, returns false
    """
    def overlapIntergene(self,gene):
        accession,start,end,strand = gene
        ints = intervals.Intervals()
        ints.setIntervals(self.intergeneDict[(accession,strand)])
        region = gene[1],gene[2]
        return region in ints
        
if __name__=="__main__":
    import unittest
    class TestIntervals(unittest.TestCase):
        def setUp(self):
            input_reads = '>NC_014561-ign-0:4778-4951+ NC_014561 4778-4951 +\n'\
                          'CGGCTCAGTTAATACCCTGAAATGTATTTCTGTGATAACCGGCCGTACAGTCACGTTAAA\n'\
                          'AAAAATCTATTAAGCTACTTGAAAGCGGACAACCGATCCCTATTTTAGCGATGATCGGGA\n'\
                          'ATATCATACCGGTCAGATGAACTTTTGGATGAACCTGTTCAGGAGATTATCACC\n'\
                          '>NC_014561-ign-0:5778-5951+ NC_014561 5778-5951 +\n'\
                          'CGGCTCAGTTAATACCCTGAAATGTATTTCTGTGATAACCGGCCGTACAGTCACGTTAAA\n'\
                          'AAAAATCTATTAAGCTACTTGAAAGCGGACAACCGATCCCTATTTTAGCGATGATCGGGA\n'\
                          'ATATCATACCGGTCAGATGAACTTTTGGATGAACCTGTTCAGGAGATTATCACC\n'
            self.input_file = "test.fa"
            read_handle = open(self.input_file,'w')
            read_handle.write(input_reads)
            read_handle.close()
        def cleanUp(self):
            os.remove(self.input_file)
                        
        def testDictionary(self):
            testHandler = IntergeneHandler(self.input_file)
            testHandler.getIntervals()
            self.assertEquals(len(testHandler.intergeneDict[('NC_014561',"+")]),2)
            
        def testContains(self):
            testHandler = IntergeneHandler(self.input_file)
            testHandler.getIntervals()
            self.assertTrue(testHandler.overlapIntergene(
                ('NC_014561',4700,5000,"+")
                ))
            self.assertTrue(testHandler.overlapIntergene(
                ('NC_014561',4800,4810,"+")
                ))
            self.assertTrue(testHandler.overlapIntergene(
                ('NC_014561',4700,4810,"+")
                ))
            self.assertFalse(testHandler.overlapIntergene(
                ('NC_014561',4700,4810,"-")
                ))
            self.assertFalse(testHandler.overlapIntergene(
                ('NC_014561',4700,5000,"-")
                ))
        
    unittest.main()
