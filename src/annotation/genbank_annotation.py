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
from collections import *
from bx.intervals import *
import annotation

loc_reg = re.compile("(\d+):>?(\d+)\S\((\S)\)")
annot_reg = re.compile("([A-Za-z0-9_]+)\s(\d+)\s(\d+)\s(\S)\s(\S+)")

        
"""
Builds multiple interval trees from annotated genes
from .ptt files
"""
class AnnotationTree(object):
    def __init__(self,fname):
        self.intervalDict = defaultdict( IntervalTree )
        self.fname=fname
        pass
    """ 
    Interval tree
    CAREFUL: All trees are indexed by the first 3 words of the organism name
    This was the best way for indexing 
    """
    def build(self):
        with open(self.fname,'r') as handle:
            key = ""
            for ln in handle:
                ln = ln.rstrip()
                toks = re.split("\s+",ln)
                
                if toks[1]=="proteins":continue
                if toks[0]=="Location":continue
                if '..' in toks[-1]:  
                    key = ' '.join(toks[:3])   
                else:
                    
                    description = ' '.join(toks[8:])
                    st,end = toks[0].split('..')
                    st,end = int(st),int(end)
                    self.intervalDict[key].add(st,end,description)
    """ Find overlapping annotated genes"""
    def find(self,org,st,end):
        overlaps = self.intervalDict[org].find(st,end)
        return overlaps
    
""" Build annotation database"""
def parse(organism,genbank_file,outHandle):
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
if __name__=="__main__":
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
            parse("yeast",self.test_file,open(self.out_file,'w'))
            
    class TestAnnotatons(unittest.TestCase):
        def setUp(self):
            test_input = test_genbank.yeast
            self.test_file = "test.gbk"
            self.out_file = "out.fa"
            handle = open(self.test_file,'w')
            handle.write(test_input)
            handle.close()
            parse("NC_12345",self.test_file,open(self.out_file,'w'))
        def tearDown(self):
            os.remove(self.test_file)
            #os.remove(self.out_file)
        def testIterator(self):
            annots = annotation.AnnotatedGenes(self.out_file)
            objs = [A for A in annots]
            self.assertEquals(len(objs),3)
            start,end,orgid,strand,locus,protid,seq = zip(*objs)
            start,end,orgid,strand,locus,protid,seq = list(start),list(end),list(orgid),list(strand),list(locus),list(protid),list(seq)
            self.assertTrue(0 in start)
            self.assertTrue("NC_12345" in orgid)
            self.assertTrue("+" in strand)
            self.assertTrue(206 in end)
    class TestAnnotationTree(unittest.TestCase):
        def setUp(self):
            self.queries = ['Halorhodospira halophila SL1, complete genome. - 1..2678452',
                            '2407 proteins',
                            'Location        Strand  Length  PID     Gene    Synonym Code    COG     Product',
                            '35..664 -       209     121588216       -       Hhal_0001       -       -       16S rRNA m(7)G-527 methyltransferase',
                            '654..2555       -       633     121588217       -       Hhal_0002       -       -       glucose inhibited division protein A',
                            '2616..3773      +       385     121588218       -       Hhal_0003       -       -       3-oxoacyl-[acyl-carrier-protein] synthase II',
                            '3787..4653      +       288     121588219       -       Hhal_0004       -       -       aminodeoxychorismate lyase apoprotein',
                            '4692..5675      +       327     121588220       -       Hhal_0005       -       -       aminodeoxychorismate lyase',
                            'Candidatus Liberibacter americanus str. Sao Paulo, complete genome. - 1..1195201',
                            '983 proteins',
                            'Location        Strand  Length  PID     Gene    Synonym Code    COG     Product',
                            '100..1596       +       498     557715459       dnaA    lam_001 -       -       ATPase involved in DNA replication initiation',
                            '1947..3122      -       391     557715460       hemN    lam_002 -       -       Coproporphyrinogen III oxidase',
                            '3143..3733      -       196     557715461       -       lam_003 -       -       Xanthosine triphosphate pyrophosphatase',
                            '3812..4558      -       248     557715462       rph     lam_004 -       -       RNase PH',
                            '4801..5598      -       265     557715463       -       lam_005 -       -       Creatinine amidohydrolase']
            self.ptt = "test.ptt"
            open(self.ptt,'w').write('\n'.join(self.queries))
        def tearDown(self):
            os.remove(self.ptt)
            
        def test1(self):
            tree = AnnotationTree(self.ptt)
            tree.build()
            genes = tree.find('Candidatus Liberibacter americanus',1947, 3122)
            self.assertTrue(len(genes)==1)
            self.assertEquals('Coproporphyrinogen III oxidase',genes[0])
            
            genes = tree.find('Halorhodospira halophila SL1,',4700, 5500)
            self.assertTrue(len(genes))
            self.assertEquals('aminodeoxychorismate lyase',genes[0])
            
    unittest.main()




    
