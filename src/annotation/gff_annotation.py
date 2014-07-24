
"""
This is similar to gff.py
Will probably want to merge eventually
"""

import sys
import os
import site
import argparse
import string
import numpy
import re
import subprocess

base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for directory_name in os.listdir(base_path):
    site.addsitedir(os.path.join(base_path, directory_name))
    
import Bio
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from collections import *
from bx.intervals import *

import fasta
import annotation


""" Build annotation database"""
def parse(organism,gff_file,index_obj,outHandle):
    index = 1
    with open(gff_file,'r') as handle:
        for ln in handle:
            if ln[0]=="#": continue
            ln = ln.rstrip()
            toks = ln.split('\t')
            chrom,seqtype,st,end,strand,description = toks[0],toks[2],toks[3],toks[4],toks[6],toks[8]
            if seqtype!="CDS": continue
            sequence = index_obj.fetch(chrom,int(st),int(end))
            protid = "_" #don't have one
            locus = "_" #don't have one
            fasta_str = ">%d\t%s\t%s\t%s\t%s\t%s\t%s\n%s\n"%(index,
                                                             organism,
                                                             st,end,strand,
                                                             locus,
                                                             protid,
                                                             sequence)
            outHandle.write(fasta_str)
            index+=1
            
if __name__=="__main__":
    import unittest
    class TestParseAnnotatons(unittest.TestCase):
        def setUp(self):
            gff_entries = [
            "##gff-version 3                ",
            "##sequence-region I 1 15072434 ",
            "##sequence-region II 1 15279421",
            "##sequence-region III 1 1378380",
            "##sequence-region IV 1 17493829",
            "##sequence-region V 1 20924180 ",
            "##sequence-region X 1 17718942 ",
            "##sequence-region MtDNA 1 13794"
            '\t'.join(["I","dust","low_complexity_region","494","503",".",".",".","Note=low_complexity region"]),
            '\t'.join(["I","Genefinder","CDS","535  ","883  ",".","-","1","ID=CDS:Y74C9A.gc1"]),
            '\t'.join(["I","Genefinder","CDS","2633 ","2712 ",".","-","0","ID=CDS:Y74C9A.gc1"]),
            ]
            self.seqs = [
                            ">I",
                            "gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagc",
                            "ctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcct",
                            "aagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaa",
                            "gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagc",
                            "ctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcct",
                            "aagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaa",
                            "gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagc",
                            "ctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcct",
                            "aagcctaagcctaagcctaagcctaagcctaaaaaattgagataagaaaa",
                            "cattttactttttcaaaattgttttcatgctaaattcaaaacgttttttt",
                            "tttagtgaagcttctagatatttggcgggtacctctaattttgcctgcct",
                            "gccaacctatatgctcctgtgtttaggcctaatactaagcctaagcctaa",
                            "gcctaatactaagcctaagcctaagactaagcctaatactaagcctaagc",
                            "ctaagactaagcctaagactaagcctaagactaagcctaatactaagcct",
                            "aagcctaagactaagcctaagcctaatactaagcctaagcctaagactaa",
                            "gcctaatactaagcctaagcctaagactaagcctaagactaagcctaaga",
                            "ctaagcctaatactaagcctaagcctaagactaagcctaagcctaaaaga",
                            "atatggtagctacagaaacggtagtacactcttctgaaaatacaaaaaat",
                            "ttgcaatttttatagctagggcactttttgtctgcccaaatataggcaac",
                            "caaaaataattgccaagtttttaatgatttgttgcatattgaaaaaaaca",
                            "tttttcgggttttttgaaatgaatatcgtagctacagaaacggttgtgca",
                            "ctcatctgaaagtttgtttttcttgttttcttgcactttgtgcagaattc",
                            "ttgattcttgattcttgcagaaatttgcaagaaaattcgcaagaaatttg",
                            "tattaaaaactgttcaaaatttttggaaattagtttaaaaatctcacatt",
                            "ttttttagaaaaattatttttaagaatttttcattttaggaatattgtta",
                            "tttcagaaaatagctaaatgtgatttctgtaattttgcctgccaaattcg",
                            "tgaaatgcaataaaaatctaatatccctcatcagtgcgatttccgaatca",
                            "gtatatttttacgtaatagcttctttgacatcaataagtatttgcctata",
                            "tgactttagacttgaaattggctattaatgccaatttcatgatatctagc",
                            "cactttagtataattgtttttagtttttggcaaaactattgtctaaacag",
                            "atattcgtgttttcaagaaatttttcatggtttttcttggtcttttcttg",
                            "gtatttttttgacaaaaatttttgtttcttgattcttgcaaaaatttttc",
                            "cgtttgacggccttgatgtgcactaccttcgcttaaatactacattttct",
                            "gaaaatgttataatagtgttcattgtttcatacaaatacttatttaatag",
                            "tatttctggttatataatttgtataaaaagtggttgacataacaaggctg",
                            "acgaaactttgtgatggctgaaaatattttcctagctttattgattttta",
                            "tttatacgtgtttgaataacttggccaaatcgccgagaaggaatagaata",
                            "ctggacgacattgtacatattttccaaaaaatcagaaagtagatgacggg",
                            "accaattctttctgtcaggttttacaaccgcccagtgcgtctacgtcaca",
                            "tgttgtataaatggttgtaaacaatatgcggaaacaatcaaatgcattcc",
                            "cataaggcataatatagaggctacaggcaatgagtatcgctctttgcttt",
                            "gtttaaagggggagtagagtttgtggggaaatatatgtttctgactctaa",
                            "ttttgcccctgataccgaatatcgatgtgaaaaaatttaaaaaaatttcc",
                            "ctgattttatattaatttttaaaatccgaaaatccattggatgcctatat",
                            "gtgagtttttaaacgcaaaattttcccggcagagacgccccgcccacgaa",
                            "accgtgccgcacgtgtgggtttacgagctgaatattttccttctattttt",
                            "atttgattttataccgattttcgtcgatttttctcattttttctcttttt",
                            "tttggtgttttttattgaaaattttgtgattttcgtaaatttattcctat",
                            "ttattaataaaaacaaaaacaattccattaaatatcccattttcagcgca",
                            "aaatcgactggagactaggaaaatcgtctggagatagaacggatcaacaa",
                            "gattattattatatcattaataatatttatcaattttcttctgagagtct",
                            "cattgagactcttatttacgccaagaaataaatttaacattaaaattgtt",
                            "catttttgaaaaaaaaataattaaaaaaacacattttttggaaaaaaaaa",
                            "taaataaaaaaaattgtcctcgaggatcctccggagcgcgtcgaatcaat",
                            "gtttccggaactctgaaaattaaatgtttgtatgattgtagaaccctttc",
                            "gctattgagatttgataacttttaagtaataaaattttcgcagtaagaca",
                            "ttaaaacatttcacaattaagctggttctgaactgtgtgaagtatattga",
                            "aaaaaactaactgatacaaaaatataattttatgatagttttctggatgt",
                            "cccaatataaacgatgtcaattctgcgacatgctacagtcatccacgaaa",
                            "gtaacccgaataccgacaaaagaagaggaacgccaactttggatagacgc",
                            "tctaggggctgattttggtcggaaaatagtcgggaaaaaatagaggacat",
                            "tacagatgaggatgaggatgaagatagaaatttgccgacaacttcgtcat"
                         ]
            
            self.gff = "test.gff"
            self.out_file = "test.out"
            self.fasta = "test.fa"
            self.faidx = "test.fai"
            open(self.fasta,'w').write('\n'.join(self.seqs))
            self.indexer = fasta.Indexer(self.fasta,self.faidx)
            self.indexer.index()
            self.indexer.load()
            open(self.gff,'w').write('\n'.join(gff_entries))
            parse("NC_12345",self.gff,self.indexer,open(self.out_file,'w'))
        def tearDown(self):
            os.remove(self.gff)    
            os.remove(self.fasta)
            os.remove(self.faidx)
            os.remove(self.out_file)
        def testIterator(self):
            annots = annotation.AnnotatedGenes(self.out_file)
            objs = [A for A in annots]
            self.assertEquals(len(objs),2)
            start,end,orgid,strand,locus,protid,seq = zip(*objs)
            start,end,orgid,strand,locus,protid,seq = list(start),list(end),list(orgid),list(strand),list(locus),list(protid),list(seq)
            self.assertTrue(535 in start)
            self.assertTrue("NC_12345" in orgid)
            self.assertTrue("-" in strand)
            self.assertTrue(883 in end)
        
            self.assertTrue(2633 in start)
            self.assertTrue("NC_12345" in orgid)
            self.assertTrue("-" in strand)
            self.assertTrue(2712 in end)
    unittest.main()
        
        