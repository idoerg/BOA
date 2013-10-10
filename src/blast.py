"""
A blast pipeline to blast bacteriocin proteins against intergenic regions


Word of caution: blastx is not thread safe :(
"""
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIStandalone

import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess

def addArgs(parser):
    parser.add_argument(\
        '--output-file', type=str, required=True,
        help='The output file containing the BLAST output')
    parser.add_argument(\
        '--keep-tmp', action='store_const', const=True, default=False,
        help='Keeps temporary files such as blast database and blast output xml')
    parser.add_argument(\
        '--num-threads', type=int, required=False, default=1,
        help='The number of threads to be run on blast')

class Record():
    """
    Records information for each blast entry that is important for multiple alignment
    """
    def __init__(self,
                 record_number,
                 description,
                 expected_value,
                 score,
                 query,             #Sequence from database
                 query_start,
                 query_end,
                 sbjct,             #Sequence that was blasted
                 sbjct_start,
                 sbjct_end):
        self.record_number = record_number
        self.description = description
        self.score = score
        self.expected_value = expected_value
        self.query = query
        self.query_start = query_start
        self.query_end = query_end
        self.sbjct = sbjct
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end

    def __str__(self):
        string = '****Alignment****\n' \
          +'Record%d\n  '%self.record_number \
          +'sequence%s\n: '+self.description \
          +'e value:%lf\n'%self.expected_value \
          +'Score:%lf\n'%self.score \
          +'Query:\t%s\t%s\t%s\n'%(self.query_start,self.query,self.query_end) \
          +'Sbjct:\t%s\t%s\t%s\n'%(self.sbjct_start,self.sbjct,self.sbjct_end)
        return string


class BLAST(object):
    def __init__(self,bacteriocins_file,intergene_file,intermediate):
        self.pid = os.getpid() #Use current pid to name temporary files
        self.protein_db = bacteriocins_file
        self.blastxml = "%s/%d.xml"%(intermediate,self.pid)
        self.genomic_query = intergene_file
        self.intermediate = intermediate
    def getFile(self):
        return self.blastxml
    """
    Build database using intergenic regions
    """
    def buildDatabase(self,base="nucleotide"):
        #cmd="formatdb -i %s -p T -o T"%(self.protein_db)
        char = "F" if base=="nucleotide" else "T"
        cmd="formatdb -i %s -p %s -o T"%(self.protein_db,char)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()

    """
    Blast sequences
    cmd: any blastall command (e.g. blastn, tblastn, blastx, blastp)
    """
    def run(self,blast_cmd="blastn",num_threads=1):
        outHandle = open(self.blastxml,'w')
        cmd="blastall -p %s -d %s -i %s -m 7 -o %s -a %d"%(blast_cmd,self.protein_db,self.genomic_query, self.blastxml, num_threads)
        proc = subprocess.Popen(cmd,shell=True)
        proc.wait()

    def cleanup(self):

        os.system("rm %s.*"%self.protein_db)
        os.remove(self.blastxml)

    def parseBLAST(self):
        input_file = self.blastxml
        hits = []
        handle = open(input_file,'r')
        try:
            blast_hits = NCBIXML.parse(handle)
            blast_records = list(blast_hits)
            i = 0
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect<0.04:
                            record=Record(record_number = i,
                                          description = alignment.title,
                                          expected_value = hsp.expect,
                                          score = hsp.score,
                                          query = hsp.query,
                                          query_start = hsp.query_start,
                                          query_end = hsp.query_end,
                                          sbjct = hsp.sbjct,
                                          sbjct_start = hsp.sbjct_start,
                                          sbjct_end = hsp.sbjct_end)
                            hits.append(record)
        except:
            print>>sys.stderr,"No blast hits"
        return hits

if __name__=="__main__":
    import unittest
    import test
    class TestNuc2NucBlast(unittest.TestCase):
        def setUp(self):
            self.refseq = "AGCTGGCGGCGCGAGGAAGAGGAACGTAGCTGGCGGCGCGAGGAAGAGGAACGT"
            self.fasta = "test.fa"
            self.refid = "test"
            test.createFasta(self.fasta,self.refid,self.refseq)
        def tearDown(self):
            try:
                os.remove(self.fasta)
            except:
                print "Nothing to remove"

        def test1(self):
            seq_obj = Seq(self.refseq)
            record = SeqRecord(seq_obj,id="YP_025292.1", name="HokC",description="toxic membrane protein, small")

            blast_obj = BLAST(self.fasta,self.fasta)
            blast_obj.buildDatabase("nucleotide")
            blast_obj.run(blast_cmd="blastn",num_threads=1)
            records = blast_obj.parseBLAST()
            self.assertEquals( len(records),1 )
            blast_obj.cleanup()

    class TestNuc2ProBlast(unittest.TestCase):
        def setUp(self):
            self.proseq = "MAIVMGR*KGAR*"
            self.nucseq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
            self.nucfasta = "nuctest.fa"
            self.profasta = "protest.fa"
            self.refid = "test"
            test.createFasta(self.nucfasta,self.refid,self.nucseq)
            test.createFasta(self.profasta,self.refid,self.proseq)
        def tearDown(self):
            try:
                os.remove(self.profasta)
                os.remove(self.nucfasta)
            except:
                print "Nothing to remove"
        def test1(self):
            seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
            record = SeqRecord(seq,id="YP_025292.1", name="HokC",description="toxic membrane protein, small")
            blast_obj = BLAST(self.profasta,self.nucfasta)
            blast_obj.buildDatabase("protein")
            blast_obj.run(blast_cmd="blastx",num_threads=1)
            records = blast_obj.parseBLAST()
            self.assertEquals( len(records),1 )
            blast_obj.cleanup()

    class TestPro2ProBlast(unittest.TestCase):
        def setUp(self):
            self.refseq = "MAIVMGR*KGAR*"
            self.fasta = "test.fa"
            self.refid = "test"
            test.createFasta(self.fasta,self.refid,self.refseq)
        def tearDown(self):
            try:
                os.remove(self.fasta)
            except:
                print "Nothing to remove"

        def test1(self):
            # seq = Seq("MAIVMGR*KGAR*")
            # record = SeqRecord(seq,id="YP_025292.1", name="HokC",description="toxic membrane protein, small")
            blast_obj = BLAST(self.fasta,self.fasta)
            blast_obj.buildDatabase("protein")
            blast_obj.run(blast_cmd="blastp",num_threads=1)
            records = blast_obj.parseBLAST()
            self.assertEquals( len(records),1 )
            blast_obj.cleanup()


    unittest.main()
