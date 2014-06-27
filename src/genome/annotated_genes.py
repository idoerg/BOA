"""
Given a directory, walks into all of the subdirectorys and identifies all of the annotated regions

Output format
>ORGANISM ST END STRAND DESCRIPTION
SEQUENCE
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

loc_reg = re.compile("(\d+):>?(\d+)\S\((\S)\)")
annot_reg = re.compile("([A-Za-z0-9_]+)\s(\d+)\s(\d+)\s(\S)\s(\S+)")
"""
A container to process fasta objects and obtain information for annotations
"""
class AnnotatedGenes(object):
    def __init__(self,infile):
        self.infile = infile
        self.iterator = SeqIO.parse(infile,"fasta")
    def __iter__(self):
        return self
    def next(self):
        try:
            record = next(self.iterator)
            toks = record.description.split('\t')
            #print toks
            assert len(toks)==7
            index,orgid,start,end,strand,locus,protid = toks
            start,end = int(start),int(end)
            sequence = str(record.seq)
            #Ignore the protein id. Locus tag is more general
            return start,end,orgid,strand,locus,protid,sequence
        except StopIteration as s:
            raise StopIteration()
        
    
def parseAnnotations(organism,genbank_file,outHandle):
    index = 1
    try:
        seq_record = SeqIO.parse(open(genbank_file), "genbank").next()
        for feature in seq_record.features:
            if feature.type == 'CDS':
                try: #because annotations are stupid
                    if "locus_tag" in feature.qualifiers:
                        locus = feature.qualifiers["locus_tag"][0]
                    elif "gene" in feature.qualifiers:
                        locus = feature.qualifiers["gene"][0]
                    else:
                        continue
                    protid,sequence = "_","_"
                    if "protein_id" in feature.qualifiers:
                        protid = feature.qualifiers["protein_id"][0]
                    if "translation" in feature.qualifiers:
                        sequence  = feature.qualifiers["translation"][0]
                    st,end,strand = loc_reg.findall(str(feature.location))[0]
                    fasta_str = ">%d\t%s\t%s\t%s\t%s\t%s\t%s\n%s\n"%(index,
                                                                     organism,
                                                                     st,end,strand,
                                                                     locus,
                                                                     protid,
                                                                     sequence)
                    outHandle.write(fasta_str)
                    index+=1
                except KeyError as k:                    
                    continue
                
    except Exception as e:
        print "Exception",e
        

    #SeqIO.write(sequences, outHandle, "fasta")    

def go(root_dir,output_file):
    outHandle = open(output_file,'w')
    for root, subFolders, files in os.walk(root_dir):
        for fname in files:
            genome_files = []
            organism,ext = os.path.splitext(os.path.basename(fname))
            if ext==".gbk":
                absfile=os.path.join(root,fname)
                parseAnnotations(organism,absfile,outHandle)
                
if __name__=="__main__":
    parser = argparse.ArgumentParser(description=\
        'Finds intergenic regions from genback files')
    parser.add_argument(\
        '--root-dir', type=str,required=False,default="",
        help='Root directory of all of the files of interest')
    parser.add_argument(\
        '--output-file', type=str, required=False,
        help='The output file containing the tab-delimited output')
    parser.add_argument(\
        '--test', action='store_const', const=True, default=False,
        help='Run the unittests')
    args = parser.parse_args()

    if not args.test:
        go(args.root_dir,args.output_file)
    else:
        del sys.argv[1:]
        import test_genbank
        import unittest
        class TestParseAnnotatons(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
            def tearDown(self):
                os.remove(self.test_file)
                #os.remove(self.out_file)
            def testParse(self):
                parseAnnotations("yeast",self.test_file,open(self.out_file,'w'))
                
        class TestAnnotatons(unittest.TestCase):
            def setUp(self):
                test_input = test_genbank.yeast
                self.test_file = "test.gbk"
                self.out_file = "out.fa"
                handle = open(self.test_file,'w')
                handle.write(test_input)
                handle.close()
                parseAnnotations("NC_12345",self.test_file,open(self.out_file,'w'))
            def tearDown(self):
                os.remove(self.test_file)
                #os.remove(self.out_file)
            def testIterator(self):
                annots = AnnotatedGenes(self.out_file)
                objs = [A for A in annots]
                self.assertEquals(len(objs),3)
                start,end,orgid,strand,locus,seq = zip(*objs)
                start,end,orgid,strand,locus,seq = list(start),list(end),list(orgid),list(strand),list(locus),list(seq)
                self.assertTrue(0 in start)
                self.assertTrue("NC_12345" in orgid)
                self.assertTrue("+" in strand)
                self.assertTrue(206 in end)
            
        unittest.main()
    
