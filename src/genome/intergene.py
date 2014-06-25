#!/usr/bin/env python
"""
Finds intergenic regions in a bacterial genome
"""
import sys
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
import os

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
#import intervals

import genbank
import blast
import intergene

from bx.intervals import *


base_path = os.path.dirname(os.path.abspath(__file__))
for directory_name in ['test']:
    site.addsitedir(os.path.join(base_path, directory_name))

loc_reg = re.compile("([A-Za-z0-9_]+)-ign-\d+:(\d+)-(\d+)(\S)")

class IntergeneHandler:
    def __init__(self,intergene_file):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.intergene_file = intergene_file
        
    """
    Retrieves intervals that can be used to determine if bacteriocin overlap any intergenic regions
    """
    def getIntervals(self):
        self.intergeneDict = defaultdict( IntervalTree )
        for record in SeqIO.parse(open(self.intergene_file,'r'),"fasta"):
            #print record.id
            #print loc_reg.findall(record.id)
            
            match = loc_reg.findall(record.id)[0]
            accession,start,end,strand = match
            self.intergeneDict[(accession,strand)].add(int(start),int(end))
            #self.intergeneDict[(accession,strand)].append( (int(start),int(end),strand) )            
    """
    Returns true if bacteriocin interval overlaps intergene region
    Otherwise, returns false
    """
    def overlapIntergene(self,gene):
        accession,start,end,strand = gene
        tree = self.intergeneDict[(accession,strand)]
        overlaps = tree.find(start,end)
        return len(overlaps)>0
        # #ints = self.intergeneDict[(accession,strand)]
        # ints = intervals.Intervals()
        # ints.setIntervals(self.intergeneDict[(accession,strand)])
        # region = gene[1],gene[2]
        # overlap = ints.search(region)
        # print str(gene), str(overlap)
        # return overlap!=None
        


# Copyright(C) 2009 Iddo Friedberg & Ian MC Fleming
# Released under Biopython license. http://www.biopython.org/DIST/LICENSE
# Do not remove this comment
def get_interregions(genbank_path,output_file,intergene_length=1,write_flag = 'w'):
    try:
        seq_record = SeqIO.parse(open(genbank_path), "genbank").next()
        cds_list_plus = []
        cds_list_minus = []
        intergenic_records = []
        # Loop over the genome file, get the CDS features on each of the strands
        for feature in seq_record.features:
            if feature.type == 'CDS':
                mystart = feature.location._start.position
                myend = feature.location._end.position
                if feature.strand == -1:
                    cds_list_minus.append((mystart,myend,-1))
                elif feature.strand == 1:
                    cds_list_plus.append((mystart,myend,1))
                else:
                    sys.stde1rr.write("No strand indicated %d-%d. Assuming +\n" %
                                      (mystart, myend))
                    cds_list_plus.append((mystart,myend,1))

        for i,pospair in enumerate(cds_list_plus[1:]):
            # Compare current start position to previous end position
            last_end = cds_list_plus[i][1]
            this_start = pospair[0]
            strand = pospair[2]
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                strand_string = "+"
                
                intergenic_records.append(
                      SeqRecord(intergene_seq,id="%s-ign-%d:%d-%d%s" % (seq_record.name,i,last_end+1,this_start,strand_string),
                      description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                            this_start,strand_string)))
        for i,pospair in enumerate(cds_list_minus[1:]):
            last_end = cds_list_minus[i][1]
            this_start = pospair[0]
            strand = pospair[2]
            if this_start - last_end >= intergene_length:
                intergene_seq = seq_record.seq[last_end:this_start]
                strand_string = "-"
                intergenic_records.append(
                      SeqRecord(intergene_seq,id="%s-ign-%d:%d-%d%s" % (seq_record.name,i,last_end+1,this_start,strand_string),
                      description="%s %d-%d %s" % (seq_record.name, last_end+1,
                                                            this_start,strand_string)))
        SeqIO.write(intergenic_records, open(output_file,write_flag), "fasta")
    except Exception as e:
        print "Error at",genbank_path,e


def go(root_dir,output_file):
    outHandle = open(output_file,'w')
    for root, subFolders, files in os.walk(root_dir):
        for fname in files:
            genome_files = []
            ext = os.path.splitext(os.path.basename(fname))[1]
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                get_interregions(absfile,output_file,
                                 intergene_length=1,
                                 write_flag='a')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
                                         'Finds intergenic regions from genback file')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--intergene-length', type=int, default=1,required=False,
        help='The path of the genbank file')
    parser.add_argument(\
        '--output-file', type=str, default='.',required=False,
        help='The fasta file containing the intergenic regions')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')

    args = parser.parse_args()

    if not args.test:
        go(args.root_dir,args.output_file)
    else:
        del sys.argv[1:]
        import unittest
        import test_genbank
        print "Testing ..."
        class TestIntergene(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
            def tearDown(self):
                os.remove(self.test_file)

            def test1(self):
                get_interregions(self.test_file,self.out_file)
                reg = re.compile("([A-Z]+\d+)-ign-\d+:(\d+)-(\d+)(\S)")
                records = [r for r in SeqIO.parse(open(self.out_file,'r'),"fasta")]
                ids = [r.id for r in records]
                #Each record has (accession,start,end,strand)
                descriptions = [reg.findall(ids[i])[0] for i in range(len(ids))]
                accession,start,end,strand = descriptions[0]
                self.assertEquals(int(start),207)
                self.assertEquals(int(end),686)

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
                            
            # def testDictionary(self):
            #     testHandler = IntergeneHandler(self.input_file)
            #     testHandler.getIntervals()
            #     self.assertEquals(len(testHandler.intergeneDict[('NC_014561',"+")]),2)
                
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
